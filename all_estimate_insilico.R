source('mixture_estimation_functions.R')

# Input seg.df (segment-wise CN values), cn.df (raw bin-wise CN values)
# theoretical mixing proportions and previously computed purity values
seg.total.df <- read.table('series_InSilico/Copynumber_segment_table.txt',
                          sep='\t', header=T, row.names=1, stringsAsFactors=F )
cn.total.df <- read.table('series_InSilico/Copynumber_raw_table.txt',
                         sep='\t', header=T, row.names=1, stringsAsFactors=F )
true.df <- read.table('series_InSilico/mixList_InSilico_true.txt', sep='\t', stringsAsFactors = F)
purity.df <- read.table('Purity_estimates_insilico.txt', sep='\t', header=T, stringsAsFactors=F)

# Parameters of subclonal ordering/labelling
# Using method sd and epsilon=0.05 as these work best for most datasets
filterMethod <- 'sd'
cutOffVec <- seq(0.05,0.35,by=0.005)
epsilon <- 0.05

###############
# Run many tests in a loop by randomly sampling a dataset in each turn
final.medians.total <- data.frame(matrix(vector(),ncol=12))
fitInfoList <- list()
segDfList <- list()

for (maxSigma in 1:4){

  simInd <- 1

  while (simInd < 51){

    # Randomly sample 5 in silico samples with minimum read count according to maxSigma
    # and generate dataset together with S0
    # Since QDNAseq tables are normalised, have to mulitple by 2
    seg.df <- seg.total.df[,c('S0',sample(names(seg.total.df)[1:(30*maxSigma)],5))]
    cn.df <- cn.total.df[,names(seg.df)]

    # Re-normalise segments to bring mode at CN=2
    reNorm <- centreSegs(seg.df)
    seg.df <- as.data.frame(t(t(seg.df)/reNorm)*2)
    cn.df <- as.data.frame(t(t(cn.df)/reNorm)*2)

    # Re-define segment boundaries across all samples and filter for segment length
    segchange <- sapply(1:(nrow(seg.df)-1),function(i) sum(seg.df[i,]!=seg.df[(i+1),])>0)
    seg.data <- data.frame(start=c(1,which(segchange)+1),
                         end=c(which(segchange),length(segchange)+1))
    seg.data$length <- seg.data$end - seg.data$start+1
    seg.sub <- subset(seg.data, length>120)
    # Curate segment values of filtered segments with normal distribution and update seg.df
    seg.tmp <- getNormalFitSegments(seg.sub, cn.df)
    seg.av <- seg.tmp[[1]]
    names(seg.av) <- names(seg.df)
    seg.df.upd <- data.frame(matrix(NA, ncol=ncol(seg.df), nrow=nrow(seg.df)))
    names(seg.df.upd) <- names(seg.df); row.names(seg.df.upd) <- row.names(seg.df)
    for(i in 1:nrow(seg.sub)){
        seg.df.upd[seg.sub$start[i]:seg.sub$end[i],] <- seg.av[i,]
    }

    # Define purity
    pHat.df <- data.frame(matrix(vector(),ncol=ncol(seg.df.upd), nrow=2))
    names(pHat.df) <- names(seg.df.upd)

    # since purity has been already computed, obtain values from either previous computation or theoretical values
    #pHat <- as.numeric(true.df[match(names(pHat.df), true.df$V3),'V1']) # use theoretical values
    #pHat[1] <- 0.245
    pHat <- as.numeric(purity.df[match(names(pHat.df), purity.df$time),'purity'])

    # Correct segment CN values with purity
    seg.df.corr <- t(t(seg.df.upd-2)*1/pHat)+2
    seg.av.corr <- as.data.frame(t(t(seg.av-2)*1/pHat)+2)

    # Discard samples with purity < 10% and only proceed if at least 3 samples (plus baseline) remain
    pth <- 0.1
    abovePurTh <- pHat>pth
    seg.av.corr <- seg.av.corr[,abovePurTh]
    if (is.null(dim(seg.av.corr))){
      next
    }
    if (ncol(seg.av.corr) < 4){
    	next
    }

    # Compute relative segment CN values according to baseline S0
    # Generate all possible permutations of the non-baseline samples
    # and identify best sample order and fitting subclonal segments for a range of clonal cutoff values
    baseSample <- 'S0'
    seg.rel <- seg.av.corr - seg.av.corr[,baseSample]
    segDfList[[paste0(maxSigma,':',simInd)]] <- seg.av.corr
    seg.rel.nonBase <- seg.rel %>% select(-one_of(baseSample))
    colToUse <- names(seg.rel.nonBase)
    nCol <- length(colToUse)
    seg.rel.toOrder <- seg.rel.nonBase[,colToUse]
    ordVec <- permutations(nCol,nCol,colToUse)
    fitInfo <- list()
    for (cutOff in cutOffVec){
        seg.rel.Eval <- filterSegmentRatios(seg.rel.toOrder, cutOff, filterMethod, 0)
        best <- findBestOrder(seg.rel.Eval, ordVec, epsilon, nCol, 0)
        fitInfo[[as.character(cutOff)]] <- best
    }
    fitInfo.df <- data.frame(cutOff = cutOffVec,
                             maxFit = sapply(fitInfo, function(x) x$max),
                             segsAboveCut = sapply(fitInfo, function(x) length(x$cons)),
                             segsInOrder = sapply(fitInfo, function(x) max(sapply(x$segs, length))))
    # Save information into list
    fitInfoList[[paste0(maxSigma,':',simInd)]] <- fitInfo

    # Identify optimal cutoff for clonal segments and compute relative resistant ratios
    # from the corresponding order and subclonal segments
    # Number of non-clonal segments to be retained: >=5 ; only proceed if at least 4 segments are classified as subclonal
    segAim <- 5
    co <- getCutOffAuto(fitInfo.df, segAim)
    fit <- fitInfo[[as.character(co)]]
    if (is.null(fit)){
      next
    }
    if(length(fit$segs[[1]])<4){
      next
    }
    ordInd <- 1 # use top order if there are multiple equivalent orders found
    seg.rel.toUse <- seg.rel[fit$segs[[ordInd]],fit$ord[ordInd,]]
    topSample <- names(seg.rel.toUse)[1]
    toEstimate <- setdiff(names(seg.rel.toUse), c(topSample,baseSample))
    final.ratios <- estimateRSegmentRatio(seg.rel.toUse, toEstimate, topSample, 1)

    # Predict relative ratios as median of all predicted values and compare to true values
    final.medians <- aggregate(final.ratios$value, by=list(final.ratios$variable), median)
    names(final.medians) <- c('time','relratio')
    final.medians$time <- as.character(final.medians$time)
    final.medians[nrow(final.medians)+1,] <- c(topSample, 1)
    final.medians$reltrue <- true.df[match(final.medians$time, true.df$V3),'V2'] / (true.df$V2[true.df$V3==topSample])

    # Predict the (absolute) subclonal ratio of all samples by fitting a mixture of Gaussian
    # and store in the same data frame
    final.medians$rat <- NA
    for (topSample in final.medians$time){
      final.medians <- estimateRGaussianFit(seg.rel.toUse, topSample, final.medians)
    }

    # Compile true values and other metrics into the result data frame
    final.medians$true <- true.df[match(final.medians$time, true.df$V3),'V2']
    topSample <- names(seg.rel.toUse)[1]
    final.medians$topratio <- final.medians[final.medians$time==topSample,'rat']
    final.medians$toptrue <- true.df[true.df$V3==topSample,'V2']
    final.medians$purity <- as.numeric(pHat[match(final.medians$time,names(pHat.df))])
    final.medians$puritytrue <- true.df[match(final.medians$time, true.df$V3),'V1']
    final.medians$maxSigma <- maxSigma
    final.medians$simInd <- simInd
    final.medians$segAim <- segAim
    final.medians$segNum <- nrow(seg.rel.toUse)
    final.medians.total <- rbind(final.medians.total, final.medians)

    simInd <- simInd+1
  }
}

write.table(final.medians.total, file='All_estimates_insilico_noise_total.txt', sep='\t', row.names=F, quote=F)
#save('fitInfoList','segDfList', file='All_estimates_insilico_noise_total.RData' )

