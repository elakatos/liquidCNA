library(pracma); library(ggplot2); library(ggpubr); library(reshape2); library(mixtools)
library(fitdistrplus); library(dplyr); library(QDNAseq); require(gtools); library(gridExtra)

###########################################################################
############## FUNCTIONS FOR ESTIMATION ALGORITHM #########################
###########################################################################


# Re-centre segments so that mode (highest peak) is at CN=2
# Input:
# - seg.df = segment matrix with bin-wise values in rows and samples in columns
# Returns: re-normalisation constant for each sample in seg.df
centreSegs <- function(seg.df){
reNorm <- c()
for (i in 1:ncol(seg.df)){
  x <- seg.df[,i]
  y <- density(x, adjust=0.75)
  peax <- findpeaks(y$y)
  reNorm <- c(reNorm, y$x[peax[which.max(peax[,1]),2]])
}
return(reNorm)
}



# Evaluate purity by finding the best fit for the peaks in the copynumber segmentation density plot
# Inputs:
# p = purity value to evaluate
# w = weight of peaks at ploidy 1,2,etc. - typically high ploidy is measured less accurately and should have a lower weight
# a = bandwidth adjustment value in generating the smoothed distribution of observed values
# x = observations (with no missing values)
# Parameters nu and nd may be defined to refine the peak finding algorithm
# Returns: the summed squared distance of observed and expected peaks
evalPurityDensity <- function(p,w,a,x,maxCN=8,nu=7,nd=4){
    seg.dist <- density(x[x>1.25 & x < 4.25], adjust=a)
    z <- seg.dist$x[findpeaks(seg.dist$y,nups=nu, ndowns=nd)[,2]]
    zhat2 <- sapply(1:length(w), function(i) min(((2+p*(i-2))-z)^2))
    if (min(((2+p*(maxCN-2))-z)^2) < p){
        zhat2 = zhat2+1
    }
    return(sqrt(sum(zhat2*w)))
}

# For each segment, compute updated segment values as the mean of bins filtered to be fit by a normal
# Inputs:
# - seg.sub = table of segments with start and end index that should be evaluated
# - cn.df = raw copy number (bin-wise) values
# Returns: mean and sd of the fitted normal distribution for each segment. The mean is used downstream as the segment CN
getNormalFitSegments <- function(seg.sub, cn.df){
    seg.sd <- data.frame(matrix(NA,nrow=nrow(seg.sub),ncol=ncol(cn.df)))
    seg.av <- data.frame(matrix(NA,nrow=nrow(seg.sub),ncol=ncol(cn.df)))
    for(i in 1:nrow(seg.sub)){
        for (j in 1:ncol(cn.df)){
            x <- cn.df[seg.sub$start[i]:seg.sub$end[i],j]
            x.fit <- fitdist(x, "norm")
            x.sub <- x[abs(x-x.fit$estimate['mean'])<2.5*x.fit$estimate['sd']]
            # iterate after filtering out outlier bins
            x.fit2 <- fitdist(x.sub, "norm")
            seg.av[i,j] <- x.fit2$estimate['mean']
            seg.sd[i,j] <- x.fit2$estimate['sd']
        }
    }
    return(list(seg.av=seg.av, seg.sd=seg.sd))
}

# Filter relative segment values (compared to baseline) to discard clonal (nonchanging) segments that stay relatively constant
# Input:
# - segs = normalised, purity-corrected segment-wise CN matrix, segments as rows, samples as columns
# - cutOff = threshold value to determine whether a segment is varying enough to be not-clonal
# - method = minmax or sd, evaluate segments based on their sd or minimum/maximum values in samples
# Returns: a filtered set of a segments that have values exceeding the threshold
filterSegmentRatios <- function(segs, cutOff, method, base=1){
    segsExtended <- cbind(segs,base) # add a further column, corresponding to the baseline
    if (!(method %in% c('minmax','sd'))){
        # if method is none of the supported ones, return the original
        print('Filtering criterion not recognised')
        return(segs)
    }
    if (method=='minmax'){
        # compute the minimal and maximal ratio in each segment, and return the ones that are far enough
        thresh.min <- 1-cutOff
        thresh.max <- 1+cutOff
        rat.max <- apply(segsExtended/base,1,max)
        rat.min <- apply(segsExtended/base,1,min)
        return(segs[(rat.max>thresh.max | rat.min<thresh.min),])
    }
    if (method=='sd'){
        # compute the standard deviation of each segment, and return the ones that are varied enough
        rat.sd <- apply(segsExtended,1,sd)
        return(segs[(rat.sd > cutOff),])
    }
}

# Evaluate a monotony of a sample order within a segment
# Input:
# - x = vector of segment CN values, ordered
# - th = error-threshold to allow for uncertain measurements of equal values
# - nCol = total number of samples observed and evaluated
# Returns: 1 or 0 depending if the segment is monotone ordered with threshold th
evalOrder <- function(x, th, nCol, base=1){
    diff <- x - c(x[2:length(x)],base)
    if (sum(x > base)>nCol/2)
    {orderFit <- sum(diff> -1*th)
    }else {
        orderFit <- sum(diff< th)}
    return(sum(orderFit==nCol))
    }

# For pre-filtered rows of segments, evaluate all possible sample-orders
# Input:
# - seg.ratios.Eval = normalised, purity-corrected segment matrix filtered to contain only non-clonal segments
# - ordVec = a list of orders to explore, typically an exhaustive list of all possible sample permutations
# - threshold = error threshold in evaluating monotony of a segment
# - nCol = number of samples to evaluate
# Returns: row-names of all segments considered, maximum fit gained, order with the maximum fit, segments monotone according to max fit
findBestOrder <- function(seg.ratios.Eval, ordVec, threshold, nCol, base=1){
    if (nrow(seg.ratios.Eval)==0){
        return(list(cons=row.names(seg.ratios.Eval),
                   max=NA,
                   ord=NA,
                   segs=list(NULL)))
    }
    fitTotal <- c()
    fitSegments <- list()
    for (ind in 1:nrow(ordVec)){
        ord <- ordVec[ind,]
        fitCounts <- apply(seg.ratios.Eval[,ord],1,function(z) evalOrder(z,threshold,nCol, base))
        fitTotal[[ind]] <- sum(fitCounts)/(nrow(seg.ratios.Eval))
        fitSegments[[ind]] <- row.names(seg.ratios.Eval)[fitCounts==1]
    }
    # return the order with the best fit and segments that match it
    inds <- which(fitTotal==max(unlist(fitTotal)))
    return(list(cons=row.names(seg.ratios.Eval),
                max=max(unlist(fitTotal)),
                ord=ordVec[inds,,drop=F],
                segs=fitSegments[inds]))
}

# Find the best (highest proportion of subclonal vs random segments) cutoff for clonal segments
# given a required minimal number of subclonal segments to recover (if known in advance)
# Input:
# - fitInfo.df = output of fitting at a range of clonal cutoff values
# - minNum = minimum number of subclonal (monotone) segments required to be returned
# Returns: cutoff value meeting requirement and maximising the fit
getCutOffAuto <- function(fitInfo.df, minNum){
x <- subset(fitInfo.df, segsInOrder > minNum)
if(nrow(x)>0){
    c <- x$cutOff[which.max(x$maxFit)]
}else{
    c <- fitInfo.df$cutOff[which.max(fitInfo.df$maxFit)]
}
return(c)
}

# Estimate (subclonal) ratio using the level of CNA as compared to a selected "high" sample
# Input:
# - seg.ratios = normalised, purity-corrected segment CN matrix of only subclonal segments
# - toEstimate = vector of samples to be estimated, i.e. not baseline, discarded or maximal resistant ratio samples
# - topSamples = samples with maximal resistant ratio to use as denominator in estimating relative ratios
# - w = weights of segments to control reliability, can be kept at 1 if no prior assumptions
# Returns: Data frame with estimated raw ratio values for each sample
estimateRSegmentRatio <- function(seg.ratios, toEstimate, topSamples, w){
    final.ratios <- data.frame(matrix(vector()))

    for (samp in toEstimate){
        tmp <- c()
        for (top in topSamples){
            tmp <- c(tmp,(seg.ratios[,samp])/(seg.ratios[,top]))
            }
        final.ratios <- rbind(final.ratios, data.frame(value=tmp, variable=samp, weight=w))
    
    }
    return(final.ratios)
}

###########################################################################
############## FUNCTIONS FOR SYNTHETIC DATASETS ###########################
###########################################################################

# Estimate the resistant ratio of a sample by fitting a Gaussian mixture to relative subclonal segment CNs
# Input:
# - seg.rel.toUse = normalised, purity-corrected segment CN of subclonal segments only
# - topSample = sample to be investigated - use in loop to evaluate for all samples
# - final.medians = output data frame with column "time' denoting the sample and "rat" denoting the estimated ratio and "rat_sd" its variance
# - (optional) nStates = number of DeltaCN states presented in the data, e.g. DCN={-1,1,2} -> nStates=3; this is used to compute the variance of the estimate
# Returns: updated final.medians containing values for ratio and variance of the estimate
estimateRGaussianFit <- function(seg.rel.toUse, topSample, final.medians, nStates=2){
    seg.top <- seg.rel.toUse[,topSample]
    start <- mean(seg.top[seg.top<0])
    if (is.nan(start)){
        # estimate the initial value for fitting from segments with DCN between 0 and 1
      start <- -1* mean(seg.top[seg.top>0 & seg.top<1])
  }
  if (is.nan(start)){ # if no starting point can be found, start from a default of 40%
      start <- -0.4
  }
  # By default the following DCN values are used: -2, -1, 1, 2, 3
  mixfit <- normalmixEM(as.numeric(seg.top),lambda=0.5,mean.constr = c("-2a","a","-a","-2a","-3a"), mu=c(2,1,-1,-2,-3)*start,sigma=0.1)
  # only return the obtained value if converged in less than 800 iterations
  if (length(mixfit$all.loglik)<800){
      final.medians$rat[final.medians$time==topSample] <- -1*(mixfit$mu[2])
      #sigma depends on the number of segments used and that depends on how many states they were distributed across
      final.medians$rat_sd[final.medians$time==topSample] <- (mixfit$sigma[1])/sqrt(nrow(seg.rel.toUse)/nStates)
  }
  return(final.medians)
}

# Generate synthetic CN measurements
# Input:
# - segments.df = matrix of absolute segment CNs (in rows) for resistant and sensitive tumour populations (in columns) and the length of each segment
# - params.df = parameter matrix of purity and resistant-ratio for each synthetic sample
# - randSegs = list of segmentIDs whose CN is random, i.e. not derived from purity and resistant-ratio
# - noise = magnitude of measurement noise added to bin-wise raw CN values
# Returns: matrix of measured segment CN values, mirroring QDNAseq measurements

measureSegsSynthetic <- function(segments.df, params.df, randSegs, noise){
seg.df <- data.frame(matrix(NA,ncol=nrow(params.df), nrow=sum(segments.df$length)))

for (i in 1:nrow(params.df)){
    seg.samp <- c()
    for (j in 1:nrow(segments.df)){
        # first compute the tumour-specific CN as the mixture of Ancestral/Sensitive and Resistant/Subclonal cells
        seg.cn.tum <- params.df$ratio[i]/100*segments.df$cnResistant[j] + (1-params.df$ratio[i]/100)*segments.df$cnSensitive[j]
        if (segments.df$id[j] %in% randSegs){
            # if the segment is unstable, sample CN values randomly
            seg.cn.tum <- runif(1, min=1, max=3.5)
        }
        seg.cn <- (1-params.df$purity[i]/100)*2 + params.df$purity[i]/100*seg.cn.tum
        # generate bin-wise raw CN measurements that match the length of the segment and contain noise controlled by sigma
        cn.tmp <- seg.cn + rnorm(segments.df$length[j], 0, (0.2*noise+0.015*noise*seg.cn))
        # true segment value is then computed as their mean
        seg.tmp <- rep(mean(cn.tmp),segments.df$length[j])
        seg.samp <- c(seg.samp, seg.tmp)
    }
    seg.df[,i] <- seg.samp
}
return(seg.df)
}

# Purity estimation of synthetic samples (a wrapper around evalPurityDensity)
# Input:
# - seg.df = segment CN matrix with segments in rows and samples in columns
# - params.df = true parameter values used in the generation of samples
# - simInd = index of dataset for looped data generation and analysis
# Returns: data frame with purity of each sample together with true purity and resistant-ratio and index
getPurSynthetic <- function(seg.df, params.df, simInd){
    # Parameters: weight of each CN state, smoothing kernel adjustment values and purity values to evaluate
    w = c(0.8,1,1,0.25, 0.05)
    adjVec = c(0.6,0.7,0.8,0.9,1,1.1,1.2,1.3)
    pVec = seq(0.03, 0.5, by=0.0025)

    pHat.df <- data.frame(matrix(vector(),ncol=6, nrow=5))

    # Run purity estimation for all samples and all bandwidths
    for(i in 1:5){
        x <- na.omit(seg.df[,i])
        pFits <- as.data.frame(sapply(adjVec,
         function(a) sapply(pVec,
          function(p) evalPurityDensity(p,w,a,x))))
        names(pFits) <- adjVec
        pFits$p <- pVec
        # Provide estimate which minimises mean and median
        mins <- c(pVec[which.min(apply(pFits[,1:length(adjVec)],1,mean))],
          pVec[which.min(apply(pFits[,1:length(adjVec)],1,median))])
        # Combine estimates with known theoretical values
        pHat.df[i,] <- c(mins, params.df$purity[i]/100, params.df$ratio[i]/100, i, simInd)
    }
    return(pHat.df)
}

# Derive the optimal order and subclonal set of segments on synthetic set (a wrapper around filterSegmentRatios and findBestOrder)
# Input:
# - seg.rel.nonBase = normalised, purity-corrected segment CN matrix of non-base samples
# Returns: information on the fit for all order and dataframe summarising these fits
getOrderSynthetic <- function(seg.rel.nonBase){
    colToUse <- names(seg.rel.nonBase)[1:ncol(seg.rel.nonBase)]
    nCol <- length(colToUse)
    seg.rel.toOrder <- seg.rel.nonBase[,colToUse]
    ordVec <- permutations(nCol,nCol,colToUse)

    # set the method by default to sd and epsilon to 0.05, as optimal for synthetic data
    filterMethod <- 'sd'
    cutOffVec <- seq(0.1,0.35,by=0.005)
    epsilon <- 0.05

    # for each cutoff value, evaluate clonal/subclonal/unstable segments and best order
    fitInfo <- list()
    for (cutOff in cutOffVec){
        seg.rel.Eval <- filterSegmentRatios(seg.rel.toOrder, cutOff, filterMethod, 0)
        best <- findBestOrder(seg.rel.Eval, ordVec, epsilon, nCol, 0)
        fitInfo[[as.character(cutOff)]] <- best
    }

    # arrange fit information into dataframe for downstream processing
    fitInfo.df <- data.frame(cutOff = cutOffVec,
                            maxFit = sapply(fitInfo, function(x) x$max),
                            segsAboveCut = sapply(fitInfo, function(x) length(x$cons)),
                            segsInOrder = sapply(fitInfo, function(x) max(sapply(x$segs, length))))
return(list(raw=fitInfo, df=fitInfo.df))
}

