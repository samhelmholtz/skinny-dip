getModalTriangleTripleFromDipOutput <- function(dipOutput){
    ## Example input:
    ##     dip() in C: n = 5; starting with  2N*dip = 1.n
    ## 'dip': LOOP-BEGIN: 2n*D= 1         [low,high] = [  1,  5]; l_lcm/gcm = ( 3, 3)
    ##   while(gcm[ix] != lcm[iv]) : ..
    ##   calculating dip .. (dip_l, dip_u) = (1.10204, 1) -> new larger dip 1.10204 (j_best = 2) lcm-centred triple (1,2,3)
    ## 'dip': LOOP-BEGIN: 2n*D= 1.102     [low,high] = [  3,  4]; l_lcm/gcm = ( 2, 2)
    ##   ** (l_lcm,l_gcm) = (2,2) ==> d := 1
    ## ... for which we would output the final triple (1,2,3) above

    triple <- c()
    ## Iterate over the dip output and find
    for(i in 1:length(dipOutput)){
        ## Search for a match of the dip point in the current string
        m <- gregexpr("triple \\(\\d+,\\d+,\\d+\\)", dipOutput[i], perl=TRUE)
        regexMatchResult <- regmatches(dipOutput[i], m)
        regexMatchResult <- regexMatchResult[[1]]
        if(length(regexMatchResult) > 0){
            ## Remove the '() '
            regexMatchSplitResult <- strsplit(regexMatchResult, split=",")
            stringIntegers <- gsub("\\(", "", regexMatchSplitResult[[1]])
            stringIntegers <- gsub("\\)", "", stringIntegers)
            stringIntegers <- gsub(" ", "", stringIntegers)
            stringIntegers <- gsub("triple", "", stringIntegers)
            triple <- as.integer(stringIntegers)
        }        
    }

    if(length(triple) != 3){
        stop("no triple found")
    }
    
    return(triple)
}

plotDipGraphics <- function(x){
    x <- x[order(x)]
    F <- ecdf(x)
    dt <- dip(x, full.result=T)
    dipPValue <- dip.test(x)
    d <- as.numeric(dt$dip)
    xu <- dt$xu
    xl <- dt$xl

    eps <- 0.000000000001
    xLeftOfStep <- x-eps

    ## Calculate the Fn+d and Fn-d points
    stepXPoints <- c(x,xLeftOfStep)
    stepXPoints <- stepXPoints[order(stepXPoints)]
    stepYPoints <- F(stepXPoints)
    stepYPointsPlusD <- stepYPoints+d
    stepYPointsMinusD <- stepYPoints-d

    ## Calculate the GCM of F+d in (-inf,xL]
    unimodalGcmXPoints <- unique(x[x<=xl])
    stepYPointsInUnimodalGcm <- F(unimodalGcmXPoints-eps)+d
    G <- gcmlcm(unimodalGcmXPoints,stepYPointsInUnimodalGcm, type="g")

    ## Calculate the LCM of F-d in [xU,inf)
    unimodalLcmXPoints <- unique(x[x>=xu])
    print(xu)
    stepYPointsInUnimodalLcm <- F(unimodalLcmXPoints)-d
    L <- gcmlcm(unimodalLcmXPoints, stepYPointsInUnimodalLcm, type="l")

    ## plot(stepXPoints, stepYPoints+d, type="l", lty=3, main=sprintf("GCM (red) and LCM (blue), dip=%f,p-value=%f",round(d,4),dipPValue$p.value),ylim=c(-0.25,1.25), ylab="Cumulative probability (Fn+d, Fn-d, and unimodal fit)", xlab="x")
    print(dipPValue$p.value)
    plot(stepXPoints, stepYPoints+d, type="l", lty=2, main="",ylim=c(-0.07,1.07), xlim=c(0,1.05),ylab="", xlab="",xaxt="n",yaxt="n")

    chicaneImg <- readJPEG("/home/sam/Documents/phd/projects/skinny-dip/tex/graphics/chicane.jpg")
    rasterImage(chicaneImg, 0.01, 0.65, 0.55, 1.02,interpolate=FALSE)
    
    lines(stepXPoints, stepYPoints-d, lty=2)
    lines(stepXPoints, stepYPoints, col=rgb(0,0,0,alpha=0.4))
    rect(xl,-0.2,xu,1.2,col=rgb(0,0,0,alpha=0.2), lwd=0)
    ## lines(c(xl,xl), c(-0.25,1.25), lty=1,col=5)

    ## GCM in (-inf,xL]
    points(G$x.knots,G$y.knots)
    lines(G$x.knots,G$y.knots,col=2)

    ## LCM in [xU,inf)
    points(L$x.knots,L$y.knots)
    lines(L$x.knots,L$y.knots,col=4)

    ## Straight line in [xl,xu]
    lines(c(xl,xu),c(F(xl-eps)+d, F(xu)-d),col=1)

    ima <- readPNG("/home/sam/Documents/phd/projects/skinny-dip/tex/graphics/car.png")
    rasterImage(ima, 0.26, 0.15, 0.36, 0.25,0.12)

    text(0.4,0.1,"G Convex",cex=1.5,col=2)    
    text(0.3,1.07,"(Analogy: chicane)",cex=1.5,font=3)    
    text(0.8,0.4,"G Max slope",cex=1.5)    
    text(1.01,0.88,"G Concave",cex=1.5,col=4)    
    text(0.55,0.58,expression(2*d),cex=1.5,col=1)
    text(0.49,0.43,expression(F+d),cex=1.5,col=1,srt=50)
    text(0.5,0.21,expression(F-d),cex=1.5,col=1,srt=50)
    text(0.15,0.35,expression(F),cex=1.5,col=1)
    text(0.37,0.55,expression(G),cex=1.5,col=1)
    arrows(0.6,0.25,0.6,0.35,code=2)
    arrows(0.6,0.6,0.6,0.48,code=2)
    arrows(0.15,0.3,0.15,0.05,code=2)
    arrows(0.37,0.5,0.37,0.26,code=2)
    axis(1, at=c(0.19,0.20), labels=c("",""), cex.axis=1);
    axis(1, at=c(0.44,0.56), labels=c("",""), cex.axis=1);
    axis(1, at=c(0.68,0.72), labels=c("",""), cex.axis=1);
    axis(1, at=c(0.775,0.825), labels=c("",""), cex.axis=1);
    axis(1, at=c(0.875,0.925), labels=c("",""), cex.axis=1);
    axis(1, at=c(0.195), labels=c("A"), cex.axis=2,tck=0,font.axis=2);
    axis(1, at=c(0.5), labels=c("B"), cex.axis=2,tck=0,font.axis=2);
    axis(1, at=c(0.7), labels=c("C"), cex.axis=2,tck=0,font.axis=2);
    axis(1, at=c(0.8), labels=c("D"), cex.axis=2,tck=0,font.axis=2);
    axis(1, at=c(0.9), labels=c("E"), cex.axis=2,tck=0,font.axis=2);
}

calculatePartialDerivative <- function(data, a, i, i1,i2,i3,beta,gamma,N){
    aBeta <- a %*% beta ## This is then equal to (x(i2)-x(i1))
    aGamma <- a %*% gamma ## This is then equal to (x(i3)-x(x1))

    eta <- i2-i1-(((i3-i1)*(aBeta))/(aGamma))

    if(eta == 0){
        stop("eta zero, which should never be the case")
    }

    betaiGamma <- beta[i]*gamma
    gammaiBeta <- gamma[i]*beta
    firstFraction <- (i3-i1)/N
    secondFraction <- ((a %*% (betaiGamma-gammaiBeta)))/(aGamma^2)
    derivative <- firstFraction*secondFraction

    if(eta>0){
        return(-derivative)       
    } else{
        return(derivative)
    }        
}

normvec <- function(x) {
    return(sqrt(sum(x^2)))
}

dipPursuit <- function(data, a,debug){
    n <- ncol(data)
    d <- nrow(data)

    if(d==1){
        return(a)
    }

    aOld <- replicate(d,0)
    ## a <- c(2,-1) ## rnorm(2)
    a <- a/normvec(a) ## Unit vector to start

    if(debug){
        print("Initial unit projection vector")
        print(a)
    }

    ## precision <- 0.0001
    stepSize <- 0.05
    iteration <- 1

    ## Our termination condition is when the height stops increasing
    dipOld <- 0

    while(TRUE){
        aOld <- a

        ## Project (making a univariate sample)
        p <- a %*% data
        p <- p[1,] ## take the first row so we have a normal vector

        ## order the sample, generating our ordered set x
        sortInverseMapping <- order(p)
        x <- p[sortInverseMapping]

        ## check the uniqueness condition
        if(length(unique(x)) != n & debug){
            print("warning: x-singular projection vector a...cannot guarantee continuity/differentiability. Does the data contain duplicates?")
        }

        ## Run the dip algorithm, capturing the output which we need for touching-triangle calculations
        dipOutput <- capture.output(dipResult <- dip(x, debug=TRUE, full.result=TRUE))

        ## Get the modal triangle indices from the collection of touching triangles
        modalTriangleTriple <- getModalTriangleTripleFromDipOutput(dipOutput)
        
        ## Calculate the partial derivative for all dimensions
        i1 <- modalTriangleTriple[1]
        i2 <- modalTriangleTriple[2]
        i3 <- modalTriangleTriple[3]
        inversei1 <- sortInverseMapping[i1]
        inversei2 <- sortInverseMapping[i2]
        inversei3 <- sortInverseMapping[i3]

        beta <- data[,inversei2]-data[,inversei1]
        gamma <- data[,inversei3]-data[,inversei1]
        
        gradient <- replicate(d,0)
        for(i in 1:d){
            gradient[i] <- calculatePartialDerivative(data,a,i, i1,i2,i3, beta,gamma,n)            
        }

        if(dipResult$dip < dipOld){
            ## We've passed our maximum, abort
            break
        }

        ## Here we know we're still on the way to the top
        ## Let's update the dip reference, and update the gradients
        dipOld <- dipResult$dip
        for(i in 1:d){
            a[i] <- a[i]+(stepSize*gradient[i])           
        }
        
        ## Normalize a
        a <- a/normvec(a);
        iteration <- iteration +1
    }

    if(debug){
        print(sprintf("Solution (in %d iterations):",iteration))
        print(a)
    }
    return(a)
}

plotDipSweep <- function(data, matrixWithProjectionVectorsAsColumns, plotx){
    y <- c()
    for(i in 1:ncol(matrixWithProjectionVectorsAsColumns)){
        a <- matrixWithProjectionVectorsAsColumns[,i]
        p <- a %*% data
        p <- p[1,] ## take the first row so we have a normal vector
        sortInverseMapping <- order(p)
        x <- p[sortInverseMapping]
        thedip <- dip(x)
        y <- c(y,thedip)
    }
    plot(plotx,y,type="l",main="", ylab="Dip statistic",xlab="Projection vector angle (radians)")
    rect(2.7,0,2.95,1,col=rgb(1,0,0,0.3),border=0)
    text(2.82,0.009,"Maximum dip",cex=1,srt=90)
}


## Note that here the indexes are always passed/returned in a global sense
getClustersInInterval <- function(data,indexStart,indexEnd, prefix,testingModalIntervalOnly, debug, significanceLevel){
    
    if(debug){
        print(sprintf("%sChecking interval [%d,%d]",prefix,indexStart,indexEnd))
    }
    ## Subset the data...that is, we want to recursively look at only the data between indexStart and indexEnd and search for modes in that distribution
    dataSubset <- data[indexStart:indexEnd]

    ## Run the dip test on the data subset. Todo: here we run both the dip.test and the dip, which does the procedure twice in total, but a simpler solution doesn't seem easily achievable
    dipTestResult <- dip.test(dataSubset)
    dipResult <- dipTestResult$dipFullResult
    modalIntervalLeft <- indexStart + dipResult$lo.hi[1] - 1
    modalIntervalRight <- indexStart + dipResult$lo.hi[2] - 1

    ## Check for non-significance using our significance threshold. If the result is non-significant, we'll assume we only have one cluster (unimodal)
    if(dipTestResult$p.value > significanceLevel){
        if(testingModalIntervalOnly){
            if(debug){
                print(sprintf("%sModal Interval [%d,%d] is unimodally distributed (p-value %f)...returning it as a cluster...",prefix, indexStart,indexEnd, dipTestResult$p.value))
            }            
            return(c(indexStart,indexEnd))
        } else{
            ## Here we know we're finding the "last" cluster. For the unimodal case where the mode is indeed just a point of a small interval, the dip test has the tendency to
            ## return a very small modal interval (makes sense). This is bad in our case because it means that our core cluster is typically going to be very small in relation to
            ## the others that are found. For this reason we need a mechanism for finding out what the "full" core cluster is
            ## Our mechanism: mirror the data and run the dip on it. We're sure that it's then multimodal, and the dip should find "fully" one of the modes as a larger interval
            if(debug){
                print(sprintf("%sInterval [%d,%d] is unimodally distributed (p-value %f)...the modal interval found in was [%d,%d]. Proceeding to mirror data to extract a fatter cluster here...",prefix,indexStart,indexEnd,dipTestResult$p.value, modalIntervalLeft,modalIntervalRight))
            }

            ## Get left and right-mirrored results
            leftMirroredDataSet <- mirrorDataSet(dataSubset, TRUE)
            leftMirroredDataSetResult <- dip(leftMirroredDataSet, full.result=TRUE)
            rightMirroredDataSet <- mirrorDataSet(dataSubset, FALSE)
            rightMirroredDataSetResult <- dip(rightMirroredDataSet, full.result=TRUE)
            if(leftMirroredDataSetResult$dip > rightMirroredDataSetResult$dip){
                clusterRange <- mapIndexRangeToOrderedMirroredDataIndexRangeInOriginalOrderedData(leftMirroredDataSetResult$lo.hi[1],leftMirroredDataSetResult$lo.hi[2], modalIntervalLeft, modalIntervalRight, length(dataSubset), TRUE, indexStart)
                if(debug){
                    print(sprintf("%sModal interval on the left-mirrored data was [%d,%d]...which corresponds to a cluser (which we'll return now) in the original data of [%d,%d].",prefix,leftMirroredDataSetResult$lo.hi[1],leftMirroredDataSetResult$lo.hi[2],clusterRange[1],clusterRange[2]))
                }
                return(clusterRange)
            } else{
                clusterRange <- mapIndexRangeToOrderedMirroredDataIndexRangeInOriginalOrderedData(rightMirroredDataSetResult$lo.hi[1],rightMirroredDataSetResult$lo.hi[2], modalIntervalLeft, modalIntervalRight, length(dataSubset), FALSE, indexStart)
                if(debug){
                    print(sprintf("%sModal interval on the right-mirrored data was [%d,%d]...which corresponds to a cluser (which we'll return now) in the original data of [%d,%d].",prefix,rightMirroredDataSetResult$lo.hi[1],rightMirroredDataSetResult$lo.hi[2],clusterRange[1],clusterRange[2]))
                }
                return(clusterRange)
            }            
        }            
    }

    if(debug){
        print(sprintf("%sModal interval [%d,%d], p=%g",prefix,modalIntervalLeft,modalIntervalRight, dipTestResult$p.value))
    }
    
    ## Otherwise, expand the modal interval to see if it has more than one cluster
    modalIntervalClusters <- Recall(data, modalIntervalLeft, modalIntervalRight,paste(prefix,"----",sep=""), TRUE, debug, significanceLevel)

    ## Now we need to look at the various cases.
    ## If we only have a left interval, we just need to proceed in it...there must be at least one cluster there
    ## If we only have a right interval, we just need to proceed in it...there must be at least one cluster there
    ## If we have both, we need to consider both...there COULD be one or more on either side
    leftIntervalExists <- (indexStart < modalIntervalLeft)
    rightIntervalExists <- (indexEnd > modalIntervalRight)
    if(!leftIntervalExists & !rightIntervalExists){
        stop("We found a statistical multimodality, but the modal interval is the full interval! This should never happen!")
    }
    if(!leftIntervalExists & rightIntervalExists){
        if(debug){
            print(sprintf("%sInterval [%d,%d] is significantly MULTIMODAL. The modal interval [%d,%d] leaves no other points to the left, so we can continue to the right...",prefix,indexStart,indexEnd,modalIntervalLeft,modalIntervalRight))
        }
        rightClusters <- Recall(data,modalIntervalRight+1,indexEnd,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        return(c(modalIntervalClusters,rightClusters))
    }
    if(leftIntervalExists & !rightIntervalExists){
        if(debug){
            print(sprintf("%sInterval [%d,%d] is significantly MULTIMODAL. The modal interval [%d,%d] leaves no other points to the right, so we can continue to the left...",prefix,indexStart,indexEnd,modalIntervalLeft,modalIntervalRight))
        }
        leftClusters <- Recall(data,indexStart,modalIntervalLeft-1,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        return(c(modalIntervalClusters,leftClusters))
    }

    ## Otherwise, we have the general case of both intervals (left and right)
    ## Here we need to check both, including the closest cluster from the modal interval    
    if(length(modalIntervalClusters)>2){        
        ## More than one cluster in the modal interval, so include the closest for the test of each extreme        
        leftClusters <- Recall(data,indexStart,modalIntervalClusters[2],paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        rightClusters <- Recall(data,modalIntervalClusters[length(modalIntervalClusters)-1],indexEnd,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        return(c(modalIntervalClusters,leftClusters,rightClusters)) 
    } else{
        if(debug){
            print(sprintf("%sInterval [%d,%d] is significantly MULTIMODAL. The modal interval [%d,%d] is unimodal with intervals left and right of it. Checking these neighbouring intervals with the modal interval...",prefix,indexStart,indexEnd,modalIntervalLeft,modalIntervalRight))
        }
        ## Single cluster in modal interval. We hence know that there exists cluster(s) outside the modal interval. Find (them) by just focusing on the extreme intervals
        leftClusters <- Recall(data, indexStart, modalIntervalRight,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)    
        rightClusters <- Recall(data,modalIntervalLeft,indexEnd,paste(prefix,"----",sep=""), FALSE, debug, significanceLevel)
        return(c(modalIntervalClusters,leftClusters,rightClusters)) 
    }
}


mergeIntervals <- function(orderedData,intervals,debug){
    ## We first need to find the merged clusters (any overlaps are merged), such that we only have mutually-exclusive clusters
    clusterStartingIndices <- intervals[seq(1,length(intervals),2)]
    clusterEndingIndices <- intervals[seq(2,length(intervals),2)]
    clusterStartingIndicesOrderMappings <- order(clusterStartingIndices)
    clusterStartingIndicesOrdered <- clusterStartingIndices[clusterStartingIndicesOrderMappings]
    clusterEndingIndicesOrdered <- clusterEndingIndices[clusterStartingIndicesOrderMappings]
    
    clusters <- c(clusterStartingIndicesOrdered[1],clusterEndingIndicesOrdered[1])
    for(i in tail(seq_along(clusterStartingIndicesOrdered), length(clusterStartingIndicesOrdered)-1)){
        ## Get the current interval
        startInQuestion <- clusters[length(clusters)-1]
        endInQuestion <- clusters[length(clusters)]

        if(endInQuestion < clusterStartingIndicesOrdered[i]){
            clusters <- c(clusters,clusterStartingIndicesOrdered[i],clusterEndingIndicesOrdered[i])
        } else if(endInQuestion < clusterEndingIndicesOrdered[i]){
            clusters[length(clusters)] <- clusterEndingIndicesOrdered[i]
        }
    }

    ## Now we do our "consolidation" step
    ## The idea is that we merge in any "tails", "fringes" or "artefacts" that aren't
    ## statistically-significant enough to cause multimodality.
    ## How? We know that our clusters are ordered.
    ## We iterate though our clusters and perform the dip test on the range defined by successfive pairs
    ## If a pair has a non-significant multimodality, we call the entire range defined by that successive pair a single cluster
    ## print("PreliminaryClusters")
    ## print(clusters)
    consolidatedClusters <- consolidateClusters(orderedData,clusters,1,debug)
    return(consolidatedClusters)
}

consolidateClusters <- function(orderedData, clusters,index,debug){
    ## If index > length-1 done
    ## do dip
    ## If significant
    ## increment index and recurse with index++
    ## else
    ## merge and recurse with index

    numClustersLeft <- length(clusters)/2
    if(index > (numClustersLeft-1)){
        return(clusters)
    }
    startingIndex <- (index*2)-1
    endingIndex <- (index+1)*2
    startingPointIndex <- clusters[startingIndex]
    endingPointIndex <- clusters[endingIndex]
    dipTestResult <- dip.test(orderedData[startingPointIndex:endingPointIndex])
    if(dipTestResult$p.value < 0.05){
        ## significant multimodality...continue with the next index
        if(debug){
            print(sprintf("Range %d to %d is significant...we're happy with that cluster!",startingPointIndex,endingPointIndex))
        }
        return(Recall(orderedData,clusters,index+1, debug))
    } else{
        if(debug){
            print(sprintf("Range %d to %d is not significant: merging",startingPointIndex,endingPointIndex))
        }
        ## not significant...merge and repeat with the same index        
        clusters <- c(clusters[1:startingIndex],clusters[endingIndex:length(clusters)])
        return(Recall(orderedData,clusters,index, debug))        
    }      
}

## The output of this method is then a data frame, where each row corresponds to a cluster, with the objects in that cluster given by the cluster start and end indices
findCoreClusters <- function(data,significanceLevel,debug){
    if(is.unsorted(data)){
        stop("data must be sorted")
    }
    
    clustersRaw <- getClustersInInterval(data,1,length(data),"----",TRUE,debug,significanceLevel) ## Note we're saying we're testing the modal interval only here, which implies that if the distribution
    ## is NOT multi-modal, we return the full data set as the cluster. This is important, otherwise the "cluster" would just be the modal interval, which even with mirroring is smaller than the full data set. Basically, we want to say: if there's no multimodal behavior in this data, we'll just call it all a single cluster

    clusters <- mergeIntervals(data,clustersRaw, debug)
    
    clusterStarts <- c()
    clusterEnds <- c()
    clusterSizes <- c()
    clusterCentres <- c()
    for(i in seq(1,length(clusters),2)){
        clusterStarts <- c(clusterStarts,clusters[i])
        clusterEnds <- c(clusterEnds,clusters[i+1])
        clusterSizes <- c(clusterSizes,clusters[i+1]-clusters[i]+1)
        clusterCentres <- c(clusterCentres,mean(data[clusters[i]:clusters[i+1]]))
    }
    
    return(data.frame(clusterStarts=clusterStarts,clusterEnds=clusterEnds,clusterSizes=clusterSizes,clusterMeans=clusterCentres))
}

mutualInformation <- function(x1,x2){
    n <- length(x1)
    ## Freedman-Diaconis rule for histogram bin width
    bw1 <- (2 * IQR(x1)) / n^(1/3)
    bw2 <- (2 * IQR(x2)) / n^(1/3)
    ## print(sprintf("Using FD bin-widths: %f, %f",bw1,bw2))
    numBins1 <- floor((max(x1)-min(x1))/bw1)
    numBins2 <- floor((max(x2)-min(x2))/bw2)
    ## print(sprintf("Using FD bin-numbers: %d, %d",numBins1,numBins2))
    dr <- discretize2d(x1,x2,numBins1,numBins2)
    return(mi.plugin(dr))
}

generateRegularSparseGrid <- function(d,l){
    filename <- sprintf("%s/%s",tempdir(),"sparse_grid_temp.csv")
    system(sprintf("cd sgpp/examples;python python_simple.py %d %d \"%s\";",d,l,filename))
    sparseGridDataFrame <- read.csv(filename,header=FALSE)
    return(t(as.matrix(sparseGridDataFrame))) ## Transpose so that the "observations" are column vectors
}

generateSimpleRandomSampleGrid <- function(d,numPoints){
    return(matrix(runif(numPoints*d),nrow=d,ncol=numPoints))
}

generateSparseGridSizesMatrix <- function(pointsThreshold,maxDimensions){
    sizes <- matrix(replicate(30*maxDimensions, Inf),maxDimensions,30)
    for(i in 1:maxDimensions){
        for(j in 1:30){
            sg <- generateRegularSparseGrid(i,j)
            sizes[i,j] <- ncol(sg)
            if(ncol(sg)>pointsThreshold){
                break
            }
        }
    }
    return(sizes)
}

findDipPursuitStartingPoint <- function(data, existingSubspace,existingSubspaceProjections,otherClusteringProjections, numSparseGridPoints, debug){
    ## Here 'otherClusteringProjections' is a matrix where the columns are the univariate data of the projections that we already have
    
    bestvec <- c()
    bestresult <- -1
    n <- ncol(data)
    d <- nrow(data)

    newClustering=is.null(existingSubspace)

    ## We need to search over a sparse grid that is one less than the null space dimension
    if(newClustering){
        nullSpaceBasis <- diag(d)
    } else{
        nullSpaceBasis <- Null(existingSubspace)
    }
    
    nullSpaceDimension <- ncol(nullSpaceBasis)

    ## If we only have one dimension in the null subspace, then that's our vector!
    if(nullSpaceDimension == 1){
        return(as.numeric(nullSpaceBasis[,1]))
    }
    
    sparseGridDimension <- nullSpaceDimension-1
    sparseGridLevel <- findSparseGridLevelBasedOnPointsThreshold(sparseGridDimension, numSparseGridPoints)
    if(debug){
        print(sparseGridDimension)
        print(sparseGridLevel)
    }
    sg <- generateRegularSparseGrid(sparseGridDimension, sparseGridLevel)
    numSparseGridPoints <- ncol(sg)

    ## TODO: This test overwrites our sparse grid with a simple-randomly-sampled grid. Remove this to return to the use of sparse grids    
    ## sg <- generateSimpleRandomSampleGrid(sparseGridDimension, numSparseGridPoints)
    
    sgAngles <- sg*pi
    hyperSpherePoints <- getNPlus1DimensionalCartesianUnitVectorsFromNDimensionalAngles(sgAngles)
    if(debug){
        print(sprintf("Using regular sparse grid of dimension %d, level %d and node count %d",sparseGridDimension, sparseGridLevel, numSparseGridPoints))
    }
    for(i in 1:numSparseGridPoints){
        if((i%%100) == 0 & debug){
            print(sprintf("Progress: iteration %d", i))
        }
        
        a <- hyperSpherePoints[,i]
        ## Rotate the vector into the null space
        if(!newClustering){
            a <- as.numeric(nullSpaceBasis %*% a)
        }
        ## Project the data
        p <- a %*% data
        p <- p[1,] ## take the first row so we have a normal R vector

        ## calculate our objective function
        if(newClustering){
            ## For a new clustering, we choose the next projection based on the dip value and also the maximum mutual information to those already found
            dipValue <- dip(p)
            ## and calculate the maximal amount of mutual information that we have to an existing projection
            ## We take the maximum because, e.g. if we have 10 clusters and this candidate is higly independent of 9 of them, but very similar to the 10th, then that's all that matters
            maxMutualInfo <- -1
            if(!is.null(otherClusteringProjections)){
                for(j in 1:ncol(otherClusteringProjections)){
                    existingProjection <- otherClusteringProjections[,j]
                    mi <- mutualInformation(existingProjection, p)
                    if(mi > maxMutualInfo){
                        maxMutualInfo <- mi
                    }
                }           
            }
            if(maxMutualInfo<0){
                objectiveMeasure <- dipValue
            } else{
                objectiveMeasure <- dipValue/maxMutualInfo
            }
        } else{
            ## Here we are building on an existing clustering. We search for the projection which gives us with the maximum minimum mutual information to those projections already found
            minMutualInfo <- Inf
            for(j in 1:ncol(existingSubspaceProjections)){
                existingProjection <- existingSubspaceProjections[,j]
                mi <- mutualInformation(existingProjection, p)
                if(mi < minMutualInfo){
                    minMutualInfo <- mi
                }
            }
            objectiveMeasure <- minMutualInfo
        }      
        
        if(objectiveMeasure > bestresult){
            bestvec <- a
            bestresult <- objectiveMeasure
        }
    }

    if(debug){
        print(bestvec)
    }

    return(bestvec)
}


## Returns a matrix which maps dimensions and grid levels to grid sizes (number of points). Limited to d=100 and l=30.
getSparseGridSizes <- function(){
    sgsizes <- read.csv("sparseGridSizesMatrix.csv", header=FALSE)
    return(matrix(as.numeric(sgsizes),100,30))    
}

findSparseGridLevelBasedOnPointsThreshold <- function(d, pointsThreshold){
    if(pointsThreshold > 100000){
        stop("Sorry, only available for up to 100000 sparse grid points")
    }
    if(d > 100){
        stop("Sorry, only available for up to 100 dimensions")
    }

    sgsizes <- getSparseGridSizes();
    sizeRow <- sgsizes[d,]
    size <- sizeRow[sizeRow > pointsThreshold][1]
    return(match(size,sizeRow))
}

getNPlus1DimensionalCartesianUnitVectorsFromNDimensionalAngles <- function(angles){
    ## Assume angles is a matrix where each column vector is of length n, representing a direction in n+1-dimensional space
    n <- nrow(angles)
    numVectors <- ncol(angles)
    nPlus1 <- n+1
    result <- matrix(replicate(nPlus1*numVectors,1), nPlus1, numVectors)
    ## Update row-wise
    for(i in 1:nPlus1){
        ## Here we want to compute the ith dimension for each vector. It's the product of a bunch of circular functions...
        if(i > 1){
            for(j in 1:(i-1)){
                result[i,] = result[i,] * sin(angles[j,])
            }
        }
        if(i < nPlus1){
            result[i,] = result[i,]*cos(angles[i,])
        }        
    }
    return(result)
}

clusterCover <- function(listOfListOfSets){
    numCollections <- length(listOfListOfSets)
    indexOfSmallestCollection <- 1
    smallestCollection <- listOfListOfSets[[indexOfSmallestCollection]]
    smallestCollectionCount <- sum(as.numeric(lapply(smallestCollection,length)))    
    if(numCollections > 1){
        for(i in 2:numCollections){
            collectionCount <- sum(as.numeric(lapply(listOfListOfSets[[i]],length)))
            if(collectionCount < smallestCollectionCount){
                indexOfSmallestCollection <- i
                smallestCollection <- listOfListOfSets[[indexOfSmallestCollection]]
                smallestCollectionCount <- collectionCount                
            }
        }
    }

    otherCollections <- listOfListOfSets[setdiff(seq_along(listOfListOfSets),indexOfSmallestCollection)]

    clusters <- list()
    
    for(indexIntoSetOfSmallestCollection in seq_along(smallestCollection)){
        smallestCollectionSet <- smallestCollection[[indexIntoSetOfSmallestCollection]]
        setResult <- list(smallestCollectionSet)
        for(otherCollectionIndex in seq_along(otherCollections)){
            otherCollection <- otherCollections[[otherCollectionIndex]]
            locals <- list()
            for(indexIntoSetOfOtherCollection in seq_along(otherCollection)){
                setFromOtherCollection <- otherCollection[[indexIntoSetOfOtherCollection]]
                for(indexIntoSetResult in seq_along(setResult)){
                    setResultSet <- setResult[[indexIntoSetResult]]
                    setIntersection <- intersect(setFromOtherCollection, setResultSet)
                    if(length(setIntersection)>0){
                        locals <- c(locals,list(setIntersection))
                    }
                }
            }
            setResult <- locals
        }
        clusters <- c(clusters,setResult)
    }
    
    return(clusters)    
}

consolidateClusterCover <- function(projections, currentClusterCover, indexInFocus, debug){
    d <- ncol(projections)
    numClustersTotal <- length(currentClusterCover)
    
    if(indexInFocus >= numClustersTotal){
        ## end of list reached
        return(currentClusterCover)
    }

    clusterInFocus <- currentClusterCover[[indexInFocus]]
    clusterMerged <- clusterInFocus
    clustersNotMerged <- list()    
    for(i in (indexInFocus+1):numClustersTotal){
        ## Check the pair to see if it should be merged
        ## How do we do that? Simply iterate over the dimensions, taking the subset of the objects in each case, doing the dip, and testing to make sure it's significant. If we find ONE that's significant, we know these are two "worthy" clusters and we keep them. Otherwise, we "can't be sure with high probability" that they're multimodal in ANY dimension, so we choose to merge them.
        otherCluster <- currentClusterCover[[i]]
        forgetOther <- TRUE
        for(j in 1:d){
            focusClusterPointsInThisDimension <- projections[clusterInFocus,j]
            otherClusterPointsInThisDimension <- projections[otherCluster,j]
            dipTestResult <- dip.test(c(focusClusterPointsInThisDimension, otherClusterPointsInThisDimension))
            if(dipTestResult$p.value<0.05){
                forgetOther <- FALSE
                break
            }
        }

        if(!forgetOther){
            ## If we get here we know that the clusterInFocus is not significantly different in any dimension to the otherCluster
            ## We know that othercluster is smaller, so we just forget about it
            clustersNotMerged <- c(clustersNotMerged, list(otherCluster))
            
            ## What do we do about that? Well this result could mean two things: that either one of the clusters is just a small amount of "noise", positioned anywhere, OR that the two clusters
            ## really should belong together.
            ## Our solution for both cases: forget the smaller cluster
            ## If so, it means that these two clusters are not significant in ANY of the dimensions
            ## We could merge them, but it's possible that they're not significant because one is simply "noise" or "too small". In that case we don't wish to merge them.
            
            
            ## If the current cluster is the larger one, we just want to take that cluster alone
            
            ## ## If the current cluster in focus is the smaller one, we want to remove it from our clusters list completely, so we can just recall this method with it removed
            ## ## Otherwise, just do nothing and the smaller one won't get added at the end
            ## if(length(clusterMerged)<length(otherCluster)){
            ##     newClusterCover <- c(head(currentClusterCover,indexInFocus-1),currentClusterCover[(indexInFocus+1):length(currentClusterCover)])
            ##     return(Recall(projections, newClusterCover,indexInFocus, debug))
            ## }
            ## clusterMerged <- c(clusterMerged,otherCluster)
        } else{
            if(debug){
                print(sprintf("Dropping cluster (%s...), length %d because it wasn't significantly different in any dimension to cluster (%s...), length %d",paste(head(otherCluster,3),collapse=","),length(otherCluster),paste(head(clusterInFocus,3),collapse=","),length(clusterInFocus)))
            }
        }
        ## else{
        ##     ## If not, we simply add the "other one" to clustersNotMerged
        ##     clustersNotMerged <- c(clustersNotMerged, list(otherCluster))
        ## }    
        
    }

    ## Now build up the new cluster collection, up to but not including the current cluster (note this will give an empty list if indexInFocus is 1)
    ## Note also that we add the clustersNotMerged to the end...i.e. we still have to check those ones
    newClusterCover <- c(head(currentClusterCover,indexInFocus-1),list(clusterMerged),clustersNotMerged)
    return(Recall(projections, newClusterCover, indexInFocus+1,debug))
}

## Shifts to zero start, then mirrors, then returns the mirrored data
## E.g. input c(2,3,4), output c(-2,-1,0,1,2)
## Assumes an ordered input
mirrorDataSet <- function(data, mirrorLeft){
    mirroredDataSet <- data
    if(mirrorLeft){
        minValue <- min(data)
        dataShifted <- data-minValue
        dataShiftedGtZero <- dataShifted[dataShifted > 0]
        dataShiftedGtZeroMirrored <- -dataShiftedGtZero
        mirroredDataSet <- c(dataShiftedGtZeroMirrored, 0, dataShiftedGtZero)
    } else{
        maxValue <- max(data)
        dataShifted <- data-maxValue
        dataShiftedLtZero <- dataShifted[dataShifted < 0]
        dataShiftedLtZeroMirrored <- -dataShiftedLtZero
        mirroredDataSet <- c(dataShiftedLtZero, 0, dataShiftedLtZeroMirrored)
    }
    return(mirroredDataSet[order(mirroredDataSet)])
}

mapIndexRangeToOrderedMirroredDataIndexRangeInOriginalOrderedData <- function(lowerIndexToMap, upperIndexToMap, lowerFallbackIndex, upperFallbackIndex, lengthOfOriginalData, mirrorLeft, offsetIndex){
    ## Let's say our original data had a length of 2
    ## Then the mirrored data will have a length of 3
    ## The mirrored data will always have an odd length
    ## The zero point will be at lengthOfOriginalData

    if((lowerIndexToMap<lengthOfOriginalData) & (upperIndexToMap>lengthOfOriginalData)){
        return(c(lowerFallbackIndex, upperFallbackIndex))
    }    

    lowerIndexMapped <- mapIndexToOrderedMirroredDataIndexInOriginalOrderedData(lowerIndexToMap, lengthOfOriginalData, mirrorLeft)
    upperIndexMapped <- mapIndexToOrderedMirroredDataIndexInOriginalOrderedData(upperIndexToMap, lengthOfOriginalData, mirrorLeft)
    if(lowerIndexMapped > upperIndexMapped){
        return(c(upperIndexMapped, lowerIndexMapped)+offsetIndex-1)
    }
    return(c(lowerIndexMapped, upperIndexMapped)+offsetIndex-1)    
}

mapIndexToOrderedMirroredDataIndexInOriginalOrderedData <- function(indexToMap, lengthOfOriginalData, mirrorLeft){
    indexMirrored <- indexToMap-lengthOfOriginalData
    if(indexMirrored<0){
        indexMirrored <- -indexMirrored
    }
    if(mirrorLeft){
        if(indexToMap>=lengthOfOriginalData){
            return(indexToMap-lengthOfOriginalData+1)
        }else{
            return(lengthOfOriginalData-indexToMap+1)
        }
        return(indexMirrored+1)
    } else{
        if(indexToMap>lengthOfOriginalData){
            return((2*lengthOfOriginalData) - indexToMap)
        } else{
            return(indexToMap)
        }
    }
}

findOptimalSubspaceForClustering <- function(data,numSparseGridPoints,significanceLevel,debug){
    
    subspace <- NULL
    dataInOrthComplement <- data
    nullSpace <- diag(nrow(data))
    while(TRUE){

        if(ncol(nullSpace)==0){
            break
        }
        
        ## Find our next direction
        if(ncol(nullSpace)==1){
            a <- as.numeric(nullSpace[,1])
        } else{
            a <- findDipPursuitStartingPoint(dataInOrthComplement,NULL,NULL,NULL,numSparseGridPoints, debug)
        }        

        ## Gradient ascent
        if(debug){
            print(sprintf("dataInOrthComplement has %d rows",nrow(dataInOrthComplement)))
            print(sprintf("a has length %d",length(a)))
            print(a)
        }
        a <- dipPursuit(dataInOrthComplement, a, debug)
        
        ## Check its significance
        dipTestResult <- dip.test(a %*% dataInOrthComplement)
        if(dipTestResult$p.value>significanceLevel){
            break
        }

        ## Find the vector in the original data space
        ## print(a)
        ## print(nullSpace)
        vectorForSubspace <- nullSpace %*% a
        
        ## Add this vector to the subspace
        if(is.null(subspace)){
            subspace <- matrix(vectorForSubspace, nrow(data),1)
        } else{
            subspace <- cbind(subspace,matrix(vectorForSubspace, nrow(data),1))
        }

        if(debug){
            print(sprintf("Added direction to subspace. Subspace now has %d dimensions.",ncol(subspace)))
        }
        
        nullSpace <- Null(subspace)
        dataInOrthComplement <- t(nullSpace) %*% data
    }
    
    return(subspace)
}

generateSyntheticDataSet <- function(numClusters, clusteringDimensionality, noiseFraction0to1, clusterSize=200){
    numNonNoisePoints <- clusterSize*numClusters
    numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1))
    numTotalPoints <- numNonNoisePoints+numNoisePoints

    ## Build each cluster
    ## Begin with a data set of just the noise points
    ## The first "clusteringDimensionality" columns are the features, the last column is the label (in this case, all labeled zero to represent noise)
    data <- matrix(c(runif(numNoisePoints*clusteringDimensionality),replicate(numNoisePoints,0)), numNoisePoints, clusteringDimensionality+1)

    ## Calculate the number of groups of three clusters. A group is either three Gaussians or three "hyperrectangles".
    numGroups <- ceiling(numClusters/3)
    ## Uniformly distribute the group centers in [0.2,0.8] (to make sure none poke out the edge of our [0,1] noise space
    groupCenters <- matrix(runif(numGroups*clusteringDimensionality, 0.2, 0.8), numGroups, clusteringDimensionality)
    ## Standard deviation for our Gaussians
    gaussianSd <- 0.005
    ## The side half-widths for our rectangles. This is 0.25 in all but one dimension, where it's 0.005 (turns it into a "long" hyperrectangle)
    rectangleSideHalfWidths <- c(0.005,replicate(clusteringDimensionality-1,0.25)) ## rep(c(0.2,0.01),length.out=clusteringDimensionality)

    ## groupIds gives us the group to which each cluster belongs. If numGroups=2, for example and numClusters=5, then this gives us a vector c(1,1,1,2,2)    
    groupIds <- rep(1:numGroups,each=3,length.out=numClusters)
    ## groupIndexes gives us the index of each cluster in the group. If numGroups=2, for example and numClusters=5, then this gives us a vector c(1,2,3,1,2)    
    groupIndexes <- rep(c(1,2,3),length.out=numClusters)
    ## isRectangular specifies which clusters are rectangular. If numGroups=2, for example and numClusters=5, then this gives us a vector c(TRUE,TRUE,TRUE,FALSE,FALSE)    
    isRectangular <- rep(c(TRUE,FALSE),each=3,length.out=numClusters)
    for(i in 1:numClusters){

        group <- groupIds[i]
        indexInGroup <- groupIndexes[i]
        clusterCenter <- as.numeric(groupCenters[group,])
        if(isRectangular[i]){
            ## If it's a rectangular group of clusters, offset it in the second-last coordinate from the middle cluster
            if(indexInGroup==2){
                clusterCenter[1] <- clusterCenter[1]+0.03
            } else if(indexInGroup==3){
                clusterCenter[1] <- clusterCenter[1]-0.03
            }                                  
        } else{
            ## If it's a Gaussian group of clusters, offset it a coordinate from the middle cluster by six sigma
            if(indexInGroup==2){
                clusterCenter[2] <- clusterCenter[2]+6*gaussianSd
            } else if(indexInGroup==3){
                clusterCenter[2] <- clusterCenter[2]-6*gaussianSd
            }            
        }

        clusterPoints <- matrix(c(replicate(clusterSize*clusteringDimensionality,0),replicate(clusterSize,i)),clusterSize,clusteringDimensionality+1)
        noisePointsWithinBounds <- replicate(nrow(data),FALSE)
        noisePointsWithinBounds[1:numNoisePoints] <- TRUE

        if(isRectangular[i]){
            upperBounds <- clusterCenter+rectangleSideHalfWidths
            lowerBounds <- clusterCenter-rectangleSideHalfWidths
            for(j in 1:clusteringDimensionality){
                clusterPoints[,j] <- runif(clusterSize, lowerBounds[j], upperBounds[j])
                noisePointsWithinBounds <- noisePointsWithinBounds & (data[,j]>lowerBounds[j]) & (data[,j]<upperBounds[j])
            }            
        } else{
            upperBounds <- clusterCenter+gaussianSd ## Assign noise points based on 1sd bound
            lowerBounds <- clusterCenter-gaussianSd
            for(j in 1:clusteringDimensionality){
                clusterPoints[,j] <- rnorm(clusterSize, clusterCenter[j], gaussianSd)
                noisePointsWithinBounds <- noisePointsWithinBounds & (data[,j]>lowerBounds[j]) & (data[,j]<upperBounds[j])
            }
        }
        
        data[noisePointsWithinBounds,clusteringDimensionality+1] <- i ## Assign noise points falling in cluster range to that cluster
        data <- rbind(data, clusterPoints)
    }

    idsColumn <- matrix(1:nrow(data),nrow(data),1)
    data <- cbind(idsColumn,data)
    data <- data[sample(1:nrow(data),nrow(data)),] ## random permute
    return(data)
}

generateRunningExampleDataSet <- function(noiseFraction0to1=0.8){

    if(noiseFraction0to1==0.8){
        set.seed(122);
    } else{
        set.seed(195);
    }
    clusterSize <- 200;
    xVals <- c();
    yVals <- c();
    c1x <- rnorm(clusterSize, 0.9, 0.02); c1y <- rnorm(clusterSize, 0.2, 0.02); c2x <- rnorm(clusterSize, 0.8, 0.02); c2y <- rnorm(clusterSize, 0.2, 0.02); c3x <- runif(clusterSize, 0.2, 0.205); c3y <- runif(clusterSize, 0.1, 0.9); c4x <- rnorm(clusterSize, 0.7, 0.01); c4y <- rnorm(clusterSize, 0.9, 0.01); c5x <- runif(clusterSize, 0.44, 0.56); c5y <- runif(clusterSize, 0.44, 0.56); c6x <- rnorm(clusterSize, 0.8, 0.02); c6y <- rnorm(clusterSize, 0.3, 0.02);
    xVals <- c(xVals, c1x,c2x,c3x,c4x,c5x, c6x);
    yVals <- c(yVals, c1y,c2y,c3y,c4y,c5y,c6y);
    numNonNoisePoints <- length(xVals);
    numNoisePoints <- round((numNonNoisePoints*noiseFraction0to1)/(1-noiseFraction0to1));
    xVals <- c(xVals, runif(numNoisePoints));
    yVals <- c(yVals, runif(numNoisePoints));
    dataMatrix <- matrix(c(xVals, yVals), length(xVals), 2);
    
    groundTruthLabels <- replicate(length(xVals),0)
    groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.9,0.2)), length(xVals), 2, byrow=TRUE))^2))<0.04] <- 1
    groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.8,0.2)), length(xVals), 2, byrow=TRUE))^2))<0.04] <- 2
    groundTruthLabels[xVals>0.2 & xVals<0.205 & yVals>0.1 & yVals<0.9] <- 3
    groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.7,0.9)), length(xVals), 2, byrow=TRUE))^2))<0.02] <- 4
    groundTruthLabels[xVals>0.44 & xVals<0.56 & yVals>0.44 & yVals<0.56] <- 5
    groundTruthLabels[sqrt(rowSums((dataMatrix-matrix(replicate(length(xVals), c(0.8,0.3)), length(xVals), 2, byrow=TRUE))^2))<0.04] <- 6

    dataMatrix <- cbind(matrix(1:length(xVals),length(xVals),1), dataMatrix, matrix(groundTruthLabels, length(xVals), 1))
    ## random permute
    dataMatrix <- dataMatrix[sample(1:length(xVals), length(xVals)),]
    set.seed(NULL);
    return(dataMatrix);    
}

## Calculates AMI for each pair of labels.
calculateMetrics <- function(truelabelslist,predictedlabelslist){
    numCases <- length(truelabelslist)
    outputfilename <- "/home/sam/tmp/metrics.csv"
    file.remove(c("/home/sam/tmp/truelabels.csv","/home/sam/tmp/predictedlabels.csv",outputfilename))
    for(i in 1:numCases){
        write.table(matrix(truelabelslist[[i]],1),"/home/sam/tmp/truelabels.csv",append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE)
        write.table(matrix(predictedlabelslist[[i]],1),"/home/sam/tmp/predictedlabels.csv",append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE)
    }
    
    system(sprintf("/home/sam/matlab/bin/matlab -nodisplay -r \"cd '/home/sam/Documents/phd/projects/skinny-dip/experiments/scripts'; calcanmi('/home/sam/tmp/truelabels.csv','/home/sam/tmp/predictedlabels.csv','%s')\"",outputfilename),ignore.stdout=TRUE)

    res <- readLines(outputfilename)
    numResults <- length(res)
    resultsMatrix <- matrix(replicate(numResults*3,-1),numResults,3)
    for(i in 1:numResults){
        resultsMatrix[i,] <- as.numeric(unlist(strsplit(res[i],split=",")))
    }
    return(resultsMatrix)
}

## This method is our new method for finding the modes in a 1d (ordered) sample
extractModalIntervals <- function(data,significanceLevel,debug){
    if(is.unsorted(data)){
        stop("data must be sorted")
    }

    ## Find the raw clusters using our recursive approach
    ## Note: here we're saying that we're not testing a modal interval. This means that, if the full distribution is not multimodal, it will only return its estimate of the mode (not the full distribution)
    clustersRaw <- getClustersInInterval(data,1,length(data),"----",FALSE,debug,significanceLevel) 

    ## Consolidation
    clusters <- mergeIntervals(data,clustersRaw, debug)
    
    clusterStarts <- clusters[seq(1,length(clusters),2)]
    clusterEnds <- clusters[seq(2,length(clusters),2)]
    return(matrix(c(data[clusterStarts],data[clusterEnds]),nrow=length(clusterStarts),ncol=2))
}

## Recursive method for finding the bounding hypercubes for clusters in the provided data, which is assumed to be a matrix where rows are observations and columns are attributes
## Subspace is an orthonormal matrix. The columns represent the directions for our data space. It is in these directions which we take univariate projections.
## ExistingHypercube is an (m x 2) matrix (m rows, representing the dimensions of the subspace. The two entries on each row are the bounds for the edge of the hypercube in the corresponding dimension.
findClusterHypercubes <- function(subspace, existingHypercube, filteredData, significanceLevel, debug){
    subspaceDimension <- ncol(subspace)
    if(nrow(existingHypercube)>=subspaceDimension){
        ## Our hypercube is complete...return it
        return(list(existingHypercube))
    }
    if(nrow(filteredData)==0){
        ## No objects: no cluster
        return(list())
    }
    nextDimension <- nrow(existingHypercube)+1

    ## Get the next direction onto which we'll project our data
    projectionVector <- as.numeric(subspace[,nextDimension])

    ## Project the data onto that direction and sort it
    univariateProjection <- as.numeric(filteredData %*% projectionVector)
    univariateProjectionOrdered <- univariateProjection[order(univariateProjection)]

    ## Get the modal intervals along this direction. We get a matrix back where the rows are the modes, with the start/end values given in the two columns
    ## We always get at least one mode back
    modalIntervals <- extractModalIntervals(univariateProjectionOrdered, significanceLevel, debug)
    numModesFound <- nrow(modalIntervals)

    hypercubes <- list()
    for(i in 1:numModesFound){
        ## Refine our hypercube
        refinedHypercube <- rbind(existingHypercube,as.numeric(modalIntervals[i,]))
        refinedData <- filteredData[filteredData[,nextDimension]>=modalIntervals[i,1] & filteredData[,nextDimension]<=modalIntervals[i,2],]
        clusterHypercubes <- Recall(subspace, refinedHypercube, refinedData, significanceLevel, debug)
        hypercubes <- c(hypercubes,clusterHypercubes)
    }

    return(hypercubes)    
}

skinnyDipClusteringFullSpace <- function(data, significanceLevel=0.05,debug=FALSE){
    ## Hypercube-detection
    hypercubes <- findClusterHypercubes(diag(ncol(data)),matrix(nrow=0,ncol=2),data,significanceLevel,debug)

    ## Object-assignment
    labels <- replicate(nrow(data),0)
    for(i in 1:length(hypercubes)){
        hypercube <- hypercubes[[i]]
        ## Get the objects that fall within this hypercube
        objectsInHypercube <- replicate(nrow(data),TRUE)
        for(j in 1:nrow(hypercube)){
            objectsInHypercube <- (objectsInHypercube & data[,j]>=hypercube[j,1] & data[,j]<=hypercube[j,2])
        }
        labels[objectsInHypercube] <- i
    }
    return(labels)
}

getInitialTractableDataSubspace <- function(data,significanceLevel){
    ## We consider the dimensionality of the data. If it's larger than 100 dimensions, we need to do dimensionality selection.
    ## In this case, the selection involves taking only those features that have a multimodal distribution (according to the dip test)
    numDataDimensions <- ncol(data)
    if(numDataDimensions>100){
        dimensionsWithMultimodalDistributions <- c()
        for(i in 1:numDataDimensions){
            dipResult <- dip.test(as.numeric(data[,i]))
            if(dipResult$p.value < significanceLevel){
                dimensionsWithMultimodalDistributions <- c(dimensionsWithMultimodalDistributions,i)
            }
        }
        if(length(dimensionsWithMultimodalDistributions)==0){
            stop("Couldn't find any features with multimodal distributions in over 100 features. Sparse grid search is not supported for over 100 features. Aborting")
        }
        print(sprintf("Initial reduction of the dimensionality from %d to %d by selecting only those features with a multimodal distribution",numDataDimensions,length(dimensionsWithMultimodalDistributions)))
        data <- data[,dimensionsWithMultimodalDistributions,drop=FALSE]
    }
    return(data)
}

skinnyDipClusteringSubspace <- function(data, numSparseGridPoints=1000,significanceLevel=0.05,debug=FALSE){
    ## For the subspace search we need the observations as columns
    ## The subspace we get back has dimensions m*d, where m is the dimension of the original data, and d is the dimension of the subspace
    data <- getInitialTractableDataSubspace(data,significanceLevel)
    subspace <- findOptimalSubspaceForClustering(t(data),numSparseGridPoints, significanceLevel,debug)
    if(is.null(subspace)){
        print("Couldn't find any directions for a subspace that are multimodal based on our threshold. Consider trying again with more points in space (i.e. sparse grid points or similar). For now, we'll return just assume that the data set is just one cluster (i.e. a single cluster label of one)")
        return(list(NULL,replicate(nrow(data),1), data)) ## Return a single cluster (all noise)
    }
    dataInSubspace <- data %*% subspace
    labels <- skinnyDipClusteringFullSpace(dataInSubspace, significanceLevel,debug)
    print(sprintf("Found %d clusters (plus noise)", max(labels)))
    return(list(subspace,labels,dataInSubspace))
}

skinnyDipKMeansWrapper <- function(dataSetFilename, outputFilename, numSparseGridPoints, significanceThreshold, debug){
    dataFromFile <- read.csv(dataSetFilename, header=FALSE)
    dataMatrix <- as.matrix(dataFromFile[,2:(ncol(dataFromFile)-1)]) ## Remove the ID and class columns
    clusteringResult <- skinnyDipClusteringSubspace(dataMatrix, numSparseGridPoints, significanceThreshold, debug)
    predictedLabels <- clusteringResult[[2]]
    dataInSubspace <- clusteringResult[[3]]
    numNonNoiseClusters <- max(predictedLabels)
    if(numNonNoiseClusters>1){
        ## Run k-means if we have two clusters or more
        clusterCenters <- matrix(nrow=0,ncol=ncol(dataInSubspace))
        for(i in 1:numNonNoiseClusters){
            dataObjectsInCluster <- dataInSubspace[predictedLabels==i,,drop=FALSE]
            clusterCenters <- rbind(clusterCenters,colMeans(dataObjectsInCluster))
        }
        tryCatch({            
            kmeansResult <- kmeans(dataInSubspace,clusterCenters)
            predictedLabels <- kmeansResult$cluster        
        },error=function(e){
            if(grepl("empty cluster",as.character(e))){
                print("We weren't able to initialize k-means with the centers, because there were one or more empty clusters that would have resulted. Reverting to the labels found by skinny-dip")
            } else if(grepl("not distinct",as.character(e))){
                print("We weren't able to initialize k-means with the centers, because the centers weren't distinct")
            } else{
                stop(e)
            }
        })        
    }

    ## We now want to get a vector of the class labels for the points which we've assigned
    write(paste(predictedLabels,collapse=","), outputFilename)
}

sparseGridVsSimpleRandomSamplingCandidateGenerationTest <- function(optimizationFunctions, optimizationBounds, dimension, sparseGridLevel)  {
    set.seed(1)
    sg <- generateRegularSparseGrid(dimension,sparseGridLevel)
    numSparseGridPoints <- ncol(sg)
    resultFrame <- data.frame(d=c(dimension),l=c(sparseGridLevel))
    for(i in seq_along(optimizationFunctions)){
        print(optimizationFunctions[i])
        optimizationBoundsForFunction <- optimizationBounds[[i]]
        optimizationFunction <- get(optimizationFunctions[i])
        sgScaled <- (sg*(optimizationBoundsForFunction[2]-optimizationBoundsForFunction[1]))+optimizationBoundsForFunction[1]
        srs <- matrix(runif(numSparseGridPoints*dimension,optimizationBoundsForFunction[1],optimizationBoundsForFunction[2]),dimension,numSparseGridPoints)

        minSg <- Inf
        minSrs <- Inf
        for(j in 1:numSparseGridPoints){
            sgCandidate <- optimizationFunction(as.numeric(sgScaled[,j]))
            srsCandidate <- optimizationFunction(as.numeric(srs[,j]))
            if(sgCandidate < minSg){
                minSg <- sgCandidate
            }
            if(srsCandidate < minSrs){
                minSrs <- srsCandidate
            }
        }

        singleResultFrame <- data.frame(sgResult=minSg,srsResult=minSrs)
        colnames(singleResultFrame) <- c(paste(optimizationFunctions[i],"Sg",sep=""), paste(optimizationFunctions[i],"Srs",sep=""))
        resultFrame <- cbind(resultFrame, singleResultFrame)        
    }
    set.seed(NULL)
    return(resultFrame)
}
