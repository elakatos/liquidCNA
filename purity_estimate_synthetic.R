source('mixture_estimation_functions.R')


# probability of random segments
pRandom <- 0.1

# Define final data frame to store results
pHat.total <- data.frame(matrix(vector(),ncol=7))
# Repeat analysis for a range of noise (sigma) values
for (noiseLevel in c(0.5, 1, 2, 4)){
pHat.noise <- data.frame(matrix(vector(), ncol=6))

for (simInd in 1:50){

# Generate synthetic sample by sampling:
# - absolute CNs for each segment (cnResistant defined relative to cnSensitive)
# - segment lengths and a Poisson-distributed number of segments that have random CNs
# - purity and resistant-ratio parameters for each sample, with the first sample representing baseline (low resistant-ratio)
segments.df <- data.frame(id = paste0('s',1:80),
                          length = sample(12:80,80,replace=T)*10,
                          cnSensitive = sample(c(2,2,2,2,2,2,2,2,2,2,2,1,1,3,3,3,3,3,3,4,4,5), 80,replace = T))
segments.df$cnResistant <- segments.df$cnSensitive + sample(c(0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,1,1,2), 80, replace=T)
segments.df$cnResistant[segments.df$cnResistant<1] <- 1

randSegs <- sample(segments.df$id, rpois(1,lambda=80*pRandom))

params.df <- data.frame(time=paste0('t',1:5),
                        purity=c(sample(seq(4,46,by=2),5,replace = T)),
                        ratio=c(sample(c(0,0,0,0,0,1,1,1,2,3,4),1),sample(seq(5,80,by=5),4,replace=T)))

# Generate measured CN matrix by computing noise-influenced bin-wise CN values
seg.df <- measureSegsSynthetic(segments.df, params.df, randSegs, noiseLevel)
# No need to renormalise as ploidy is close to 2

# Compute purity and store results
pHat <- getPurSynthetic(seg.df,params.df,simInd)
pHat.noise <- rbind(pHat.noise, pHat)
    }

pHat.noise$noise <- noiseLevel
pHat.total <- rbind(pHat.total, pHat.noise)

}
write.table(pHat.total, file='Purity_estimates_synthetic.txt', sep='\t', row.names=F, col.names=F, quote=F)


