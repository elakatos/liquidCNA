source('mixture_estimation_functions.R')


pRandom <- 0.1

# Initialise final table to store results and lists to store intermediate information
ratHat.total <- data.frame(matrix(vector(),ncol=13))
fitInfoList <- list()
segRelList <- list()

# Iterate through 4 noise levels and make 50 simulated datasets each
for (noiseLevel in c(0.5, 1, 2, 4)){
  simInd <- 1

  while (simInd < 51){

  # Generate synthetic sample by sampling:
  # - absolute CNs for each segment (cnResistant defined relative to cnSensitive)
  # - segment lengths and a Poisson-distributed number of segments that have random CNs
  # - purity and resistant-ratio parameters for each sample, with the first sample representing baseline (low resistant-ratio)
  # Also randomise the number of samples: sampled from 2-10
    segments.df <- data.frame(id = paste0('s',1:80),
                              length = sample(12:80,80,replace=T)*10,
                              cnSensitive = sample(c(2,2,2,2,2,2,2,2,2,2,2,1,1,3,3,3,3,3,3,4,4,5), 80,replace = T))
    segments.df$cnResistant <- segments.df$cnSensitive + sample(c(0,0,0,0,0,0,0,0,0,0,0,-1,-1,-1,1,1,1,1,2,2,-2), 80, replace=T)
    segments.df$cnResistant[segments.df$cnResistant<1] <- 1

    randSegs <- sample(segments.df$id, rpois(1,lambda=80*pRandom))

    # Define or sample the number of samples to be generated and analysed
    sampleNo <- sample(2:8,1) # sample between 2 and 8
    #sampleNo <- 2 # focus on n=2 case
    print(sampleNo)
    bootstrapSingleSamp <- TRUE # also define if we want to bootstrap when n=2

    params.df <- data.frame(time=paste0('X',1:(sampleNo)),
                            purity=c(sample(seq(4,46,by=2),sampleNo,replace = T)),
                            ratio=c(sample(c(0,0,0,0,0,1,1,1,2,3,4),1),sample(seq(5,80,by=5),(sampleNo-1),replace=T)))

    # Generate measured CN matrix by computing noise-influenced bin-wise CN values
    seg.df <- measureSegsSynthetic(segments.df, params.df, randSegs, noiseLevel)
    #reNorm <- centreSegs(seg.df)
    #seg.df <- as.data.frame(t(t(seg.df)/reNorm)*2)

    # Use either estimates or known purity and correct CN profiles
    #pHat <- getPurSynthetic(seg.df,params.df,simInd)
    #pVec <- as.numeric(pHat[,1])
    pVec <- params.df$purity/100

    seg.df.corr <- t(t(seg.df-2)*1/pVec)+2
    seg.av.corr <- as.data.frame(seg.df.corr[!duplicated(seg.df.corr),])

    # Exclude cases with too low purity and discard this simulation if the base sample is not reliable or there are not enough non-base samples
    seg.av.corr <- seg.av.corr[,pVec>=0.1]
    if (is.null(dim(seg.av.corr))){
      next
    }
    if (ncol(seg.av.corr)<2){
      next
    }
    if ( !('X1' %in% names(seg.av.corr))){
      next
    }

    # Compute relative segment CN values according to baseline X1
    # Generate all possible permutations of the non-baseline samples
    # and identify best sample order and fitting subclonal segments for a range of clonal cutoff values
    baseSample <- 'X1'
    seg.rel <- seg.av.corr - seg.av.corr[,baseSample]
    #segRelList[[paste0(noiseLevel,':',simInd)]] <- seg.rel # save intermediate value for later checking if needed
    seg.rel.nonBase <- seg.rel %>% select(-one_of(baseSample))
    nCol <- ncol(seg.rel.nonBase)
    fitsAll <- getOrderSynthetic(seg.rel.nonBase)

    fitInfo <- fitsAll[[1]]
    #fitInfoList[[paste0(noiseLevel,':',simInd)]] <- fitInfo # save intermediate value for later checking if needed
    fitInfo.df <- fitsAll[[2]]

    # Identify optimal cutoff for clonal segments and compute relative resistant ratios
    # from the corresponding order and subclonal segments 
    # Number of non-clonal segments to be retained: >=12 ; only proceed if a minimum of 5 subclonal segments are found
    segAim <- 12
    co <- getCutOffAuto(fitInfo.df, segAim)
    fit <- fitInfo[[as.character(co)]]
    if (is.null(fit)){
      next
    }
    if(length(fit$segs[[1]])<5){
      next
    }
    # Here the case of n=2 (only 1 non-baseline sample) has to be handled separately if bootstrapping is to be used
    # Otherwise run standard procedure
    if (nCol==1 & bootstrapSingleSamp){
        ratiosVec <- bootstrapSegmentsAndEstimateR(seg.rel, fit, 50)
        final.medians <- data.frame(time=names(seg.rel.nonBase), relratio=1, reltrue=1,
                                    rat=mean(ratiosVec, na.rm=T), rat_sd=sd(ratiosVec, na.rm=T))
        topSample <- names(seg.rel.nonBase)[1]
    }else{
        ordInd <- 1
        seg.rel.toUse <- seg.rel[fit$segs[[ordInd]],fit$ord[ordInd,], drop=F]
        topSample <- names(seg.rel.toUse)[1]

        # Predict ratios as median of all predicted values and compare to true values
        if( nCol==1){
            final.medians <- data.frame(matrix(vector(),ncol=2))
        }else{
            toEstimate <- setdiff(names(seg.rel.toUse), c(topSample,baseSample))
            final.ratios <- estimateRSegmentRatio(seg.rel.toUse, toEstimate, topSample, 1)
            final.medians <- aggregate(final.ratios$value, by=list(final.ratios$variable), median)
        }
        names(final.medians) <- c('time','relratio')
        final.medians$time <- as.character(final.medians$time); final.medians$relratio <- as.numeric(final.medians$relratio)
        final.medians[nrow(final.medians)+1,] <- c(topSample, 1)
        final.medians$reltrue <- params.df[match(final.medians$time, params.df$time),'ratio']/params.df$ratio[params.df$time==topSample]


        # Predict the (absolute) resistant ratio of all samples by fitting a mixture of Gaussian
        # and store in the same data frame
        final.medians$rat <- NA
        final.medians$rat_sd <- NA
        tryCatch(
            for (topSample in final.medians$time){
                final.medians <- estimateRGaussianFit(seg.rel.toUse, topSample, final.medians)
            },
            error=function(e){print('Error in Gaussian fitting, continue.')}
            )
        topSample <- names(seg.rel.toUse)[1]
    }
    
    # Compile true values and other metrics into the result data frame
    final.medians$true <- params.df[match(final.medians$time, params.df$time),'ratio']/100
    final.medians$topratio <- final.medians[final.medians$time==topSample,'rat']
    final.medians$toptrue <- params.df[params.df$time==topSample,'ratio']/100
    final.medians$purity <- params.df[match(final.medians$time, params.df$time),'purity']/100
    final.medians$ratBase <- params.df[params.df$time=='X1','ratio']/100
    final.medians$simInd <- simInd
    final.medians$sigma <- noiseLevel
    final.medians$sampNum <- ncol(seg.av.corr)
    ratHat.total <- rbind(ratHat.total, final.medians)
    simInd = simInd + 1

  }
}
write.table(ratHat.total, file='Ratio_estimates_synthetic_noise_samplenumber_total.txt', sep='\t', row.names=F, quote=F)
#save('fitInfoList','segRelList', file='Ratio_estimates_synthetic_noise_samplenumber_total.RData' )


