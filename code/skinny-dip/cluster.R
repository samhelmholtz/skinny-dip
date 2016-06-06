source('load.R')
dataSetFilename <- commandArgs()[6]
outputFilename <- commandArgs()[7]

data <- read.csv(dataSetFilename,header=FALSE)
dataMatrix <- as.matrix(data[,2:(ncol(data)-1)])
startTime <- proc.time()
labels <- skinnyDipClusteringFullSpace(dataMatrix, 0.05,FALSE)
endTime <- proc.time()
runningTime <- endTime-startTime;
elapsedSeconds <- as.numeric(runningTime)[3];
write(paste(c(elapsedSeconds),collapse=","),outputFilename);
## print(sprintf("Found %d clusters (plus noise)", max(labels)))
## write(paste(labels,collapse=","), outputFilename)

## skinnyDipKMeansWrapper(dataSetFilename, outputFilename, 1000, 0.05, FALSE)
