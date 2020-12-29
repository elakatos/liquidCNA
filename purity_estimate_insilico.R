source('mixture_estimation_functions.R')

# Input seg.df (segment-wise CN values) and cn.df (raw bin-wise CN values)
seg.total.df <- read.table('series_InSilico/Copynumber_segment_table.txt',
                          sep='\t', header=T, row.names=1, stringsAsFactors=F )
cn.total.df <- read.table('series_InSilico/Copynumber_raw_table.txt',
                         sep='\t', header=T, row.names=1, stringsAsFactors=F )

# Parameters of purity estimation
w = c(0.8,1,1,0.5, 0.1, 0.025)
adjVec = c(0.6,0.7,0.8,0.9,1,1.1,1.2,1.3)
pVec = seq(0.03, 0.5, by=0.0025)

###############

# Select ALL samples to perform purity estimation
seg.df <- seg.total.df
cn.df <- cn.total.df[,names(seg.df)]

# Re-define segment boundaries across all samples and filter for segment length
segchange <- sapply(1:(nrow(seg.df)-1),function(i) sum(seg.df[i,]!=seg.df[(i+1),])>0)
seg.data <- data.frame(start=c(1,which(segchange)+1),
                     end=c(which(segchange),length(segchange)+1))
seg.data$length <- seg.data$end - seg.data$start+1
seg.sub <- subset(seg.data, length>120)

# Curate segment values of filtered segments by fitting a normal distribution and update seg.df
seg.tmp <- getNormalFitSegments(seg.sub, cn.df)
seg.av <- seg.tmp[[1]]
names(seg.av) <- names(seg.df)
seg.df.upd <- data.frame(matrix(NA, ncol=ncol(seg.df), nrow=nrow(seg.df)))
names(seg.df.upd) <- names(seg.df); row.names(seg.df.upd) <- row.names(seg.df)
for(i in 1:nrow(seg.sub)){
    seg.df.upd[seg.sub$start[i]:seg.sub$end[i],] <- seg.av[i,]
}
# Re-normalise segments to bring mode to CN=2
reNorm <- centreSegs(na.omit(seg.df.upd))
seg.df.upd <- as.data.frame(t(t(seg.df.upd)/reNorm)*2)

# Define purity
pHat.df <- data.frame(matrix(vector(),ncol=ncol(seg.df.upd), nrow=2))
names(pHat.df) <- names(seg.df.upd)
for(i in 1:ncol(seg.df.upd)){
	x <- na.omit(seg.df.upd[,i])
	pFits <- as.data.frame(sapply(adjVec,
       function(a) sapply(pVec,
                          function(p) evalPurityDensity(p,w,a,x))))
pFits$p <- pVec
mins <- c(pVec[which.min(apply(pFits[,1:length(adjVec)],1,mean))],
          pVec[which.min(apply(pFits[,1:length(adjVec)],1,median))])
pHat.df[,names(seg.df.upd)[i]] <- mins
}

output.df <- data.frame(time=names(pHat.df), purity=as.numeric(pHat.df[1,]), maxSigma=c(rep(c(1,2,3,4),each=30),1))

write.table(output.df, file='Purity_estimates_insilico.txt', sep='\t', row.names=F, quote=F)

