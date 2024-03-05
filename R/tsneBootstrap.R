tsneBootstrap <-
function(group=c("sudo","docker"), scratch.folderDOCKER, scratch.folderHOST, file, nPerm, permAtTime, percent, range1, range2, separator, logTen=0, seed=111, sp=0.8, clusterPermErr=0.05, perplexity=10){

  isDocker <- is_running_in_docker()

    if (isDocker == TRUE){
    scratch.folderHOST <- gsub("\\\\", "/", scratch.folderHOST)
    }

    if(isDocker == FALSE){
   scratch.folderDOCKER=scratch.folderHOST
   }

  permutationClustering(group=group, scratch.folderDOCKER=scratch.folderDOCKER, scratch.folderHOST=scratch.folderHOST, file=file, nPerm=nPerm, permAtTime=permAtTime, percent=percent, range1=range1, range2=range2, separator=separator, logTen=logTen, clustering="tSne", perplexity=10 , seed=seed, rK=0)
  permAnalysis(group=group, scratch.folderDOCKER=scratch.folderDOCKER, scratch.folderHOST=scratch.folderHOST, file=file, range1=range1, range2=range2, separator=separator, sp=sp, clusterPermErr=clusterPermErr, maxDeltaConfidence=0.01, minLogMean=0.05)


}
