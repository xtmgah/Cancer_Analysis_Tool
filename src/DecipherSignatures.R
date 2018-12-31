
cat('\n')
cat(paste0('Deciphering ',as.character(k_mer),'-mer mutational signatures...\n'))

#suppressMessages(library('parallel'))
#local_lib_path <- paste0(getwd(),'/../R_libs')
#.libPaths(c(.libPaths(),local_lib_path))
suppressMessages(library(doParallel))







# load the matrix M
M <- simplify2array(read.table(paste0("result/",file_name,".csv"), sep = ","))

M <- matrix(as.numeric(M),dim(M)[1])

sum_of_mutations <- colSums(M)
zero_samples <- which(sum_of_mutations == 0)
if(length(zero_samples) != 0)                 # If there are samples zero value for all mutation types, we delete them at this step.
{                                             # At the end, we again include these samples by inserting sxtra zero columns into matrix E. 
  M <- M[,-zero_samples]
}

if(k_mer == 5) {
  G <- dim(M)[2] - 1  
  M <- M[,-1]
  M <- matrix(M,ncol = G)
} else if(k_mer == 3) {
  G <- dim(M)[2]
}
len <- dim(M)[1]








#####################################################################
### Start the algorithm
#####################################################################




## Step 1 (Dimension Reduction)
# Greedy algorithm used for finding the largest subset of rows which sum up to 1% of the total sum


t <- sum(M)/100
r <- as.numeric(rowSums(M))
large <- max(r)+1
inds <- vector()
s <- 0
while(1)
{
  i <- as.numeric(which.min(r))
  if( (s + as.numeric(r[i]))  > t)
    break()
  if(length(which(colSums(M[-c(inds,i),])==0)) == 0)
  {
    s <- s + as.numeric(r[i])
    inds <- c(inds , i)
  }
  r[i] <- large
}
if(length(inds) != 0)
{
  M <- M[-inds,]  # Delete the rows...
}
M <- matrix(M, ncol = G)
remaining_mut_types <- setdiff(c(1:len),inds)  # mutation types that are present after step 1
write.table(remaining_mut_types,
            file=paste0(destination_folder,"remaining_mut_types.txt"),
            row.names = FALSE,col.names = FALSE)
rm(t,r,large,inds,s,i)
K <- dim(M)[1]

cat('\n')
cat("Dimension Reduction done\n")







# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
NMF_total_max <- 5000                                   # maximum number of NMF iterations
NMF_iters <- 1000                                       # number of NMF iterations in each epoch
NMF_max_epoches <- ceiling(NMF_total_max/NMF_iters)     # fine tune the maximum number of NMF
NMF_conv <- 1e-2                                        # stop criteria for NMF

Boot_total_max <- 100                                   # maximum number of bootstrap iterations
Boot_iters <- max(20,number_of_cpu_cores)               # number of bootstrap iterations in each epoch
Boot_max_epoches <- ceiling(Boot_total_max/Boot_iters)  # fine tune the maximum number of bootstraps
Boot_conv <- 1e-2                                       # stop criteria for bootstrap

start_N <- 1
Max_N <- 12

increment_progress_by <- 0.99/(3*Max_N+1)
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/





cat('\n')
cat(paste("With",as.numeric(number_of_cpu_cores),"CPU cores\n"))
cat(paste("NMF_iters =",as.character(NMF_iters),'\n'))
cat(paste("NMF_total_max =",as.character(NMF_total_max),'\n'))
cat(paste("NMF_conv =",as.character(NMF_conv),'\n'))
cat(paste("Boot_iters =",as.character(Boot_iters),'\n'))
cat(paste("Boot_total_max =",as.character(Boot_total_max),'\n'))
cat(paste("Boot_conv =",as.character(Boot_conv),'\n'))
cat('\n')





for(N in c(start_N:Max_N))
{
  
  nmf <- function(boot_num)
  {
    cat(paste0('\n-------------------------- Bootstrap ',as.character(boot_num),'/',as.character(Boot_iters),'\n'))
    cat('------------------------------------ Bootstrapping...')
    ## Step 2 (Bootstrap)
    bootstrap_column <- function(column)
    {
      sampled <- rep(0,K)
      samps <- sample(c(1:K),sum(M[,column]),replace=TRUE,prob=M[,column])
      ts <- table(samps)
      sampled[as.numeric(names(ts))] <- as.numeric(ts)
      return(sampled)
    }
    M.r <- sapply(c(1:G),bootstrap_column)
    cat('Done\n')
    
    
    
    ## Step 3 (NMF)
    P0 <- matrix(runif(K*N),nrow=K,ncol=N)
    E0 <- matrix(runif(N*G),nrow=N,ncol=G)
    P <- P0
    E <- E0
    U <- matrix(1,nrow=K,ncol=G)
    for(epoch in 1:NMF_max_epoches)
    {

      cat(paste0('------------------------------------ Performing NMF with ',
                 as.character(NMF_iters),' iterations : epoch ',
                 as.character(epoch),'...'))
      
      st_nmf_epoch <- Sys.time()

      max_element_change <- -1
      for(i in 1:NMF_iters)
      {
        # Normalaize matrix P and then correct matrix E so that PE remains fixed
        sums <- colSums(P)
        P <- t(t(P) / sums)
        E <- E * sums
        
        PE <- P%*%E
        zero_elements <- which(PE == 0)
        if(length(zero_elements) != 0)
        {
          mi <- min(PE[-zero_elements])
          PE[zero_elements] <- mi
        }
        E <- E * (t(P)%*%(M.r/(PE)))/(t(P)%*%U)
        
        PE <- P%*%E
        zero_elements <- which(PE == 0)
        if(length(zero_elements) != 0)
        {
          mi <- min(PE[-zero_elements])
          PE[zero_elements] <- mi
        }
        PP <- P * ((M.r/(PE))%*%t(E))/(U%*%t(E))
        
        max_element_change <- max(max_element_change,max(abs(PP-P))) 

        P <- PP
      }

      se_nmf_epoch <- Sys.time()

      cat(paste0('Done (after ',format(se_nmf_epoch-st_nmf_epoch),')\n'))

      cat(paste0('----------------------------------------------- maximum element change in this epoch = ',as.character(max(max_element_change)),'\n'))
      
      if(max(max_element_change) < NMF_conv) { break() }
    }



    # Normalaize each column of matrix P and then correct matrix E so that PE remains fixed
    sums <- colSums(P)
    P <- t(t(P) / sums)
    E <- E * sums
    
    
    return(list(P,E))
  }
  
  ## Step 4 (Iterate)
  
  incProgress(increment_progress_by, detail = paste0("\nN = ",as.character(N)," : Performing NMF"))
  #Sys.sleep(0.01)
  
  
  cat(paste0('N = ',as.character(N)),'\n')
  SP <- matrix(nrow = 0,ncol = K)
  SE <- matrix(nrow = 0,ncol = G)
  
  history_of_centroids <- list()
  
  for(epoch_boot in 1:Boot_max_epoches)
  {
    
    #--------------------------------------------------------------------
    cat(paste0('----- ','Epoch ',as.character(epoch_boot),' : With bootstrapping ',as.character(Boot_iters),' times\n'))
    
    cat('--------------- Performing NMF for each bootstrapped data\n')
    

    ########################################################
    # Save SP and SE to the current folder and delete them from environment so that it becomes light.
    if(epoch_boot != 1)
    {
      cat('--------------- pushing extra results from previous epoches...')
      st_push <- Sys.time()
      write.table(SP,'temp/SP.txt',row.names = F,col.names = F)
      write.table(SE,'temp/SE.txt',row.names = F,col.names = F)
      rm("SP")
      rm("SE")
      se_push <- Sys.time()
      cat(paste0('Done (after ',format(se_push-st_push),')\n'))
    }
    ########################################################
    

    ##### Paralelize...
    if(number_of_cpu_cores > 1)
    {
      cat('--------------- prepareing for parallelization...')
      st_prepare <- Sys.time()
      
      # cl <- makeForkCluster(number_of_cpu_cores)
      
      cl <- makeCluster(number_of_cpu_cores)
      registerDoParallel(cl)
      clusterCall(cl, function(x) .libPaths(x), .libPaths())
      
      se_prepare <- Sys.time()
      
      cat(paste0('Done (after ',format(se_prepare-st_prepare),')\n'))
      

      cat('--------------- parallelizing...')
      st_nmf <- Sys.time()
     

      # result <- parLapply(cl,c(1:Boot_iters),nmf)
      
      
      result <- foreach(i=1:Boot_iters) %dopar% {nmf(i)}
      


      stopCluster(cl)
      se_nmf <- Sys.time()

      cat(paste0('Done (after ',format(se_nmf-st_nmf),')\n'))
    } else {
      result <- lapply(c(1:Boot_iters),nmf)
    }
    #--------------------------------------------------------------------
    
    
    ########################################################
    if(epoch_boot == 1)
    {
      SP <- matrix(nrow = 0,ncol = K)
      SE <- matrix(nrow = 0,ncol = G)
    } else {
      # Retrive the SP and SE stored in the current folder
      cat('--------------- pulling the results from previous epoches...')
      st_pull <- Sys.time()
      SP <- as.matrix(read.table('temp/SP.txt'),header=F)
      SE <- as.matrix(read.table('temp/SE.txt'),header=F)
      file.remove('temp/SP.txt')
      file.remove('temp/SE.txt')
      se_pull <- Sys.time()
      cat(paste0('Done (after ',format(se_pull-st_pull),')\n'))
    }
    ########################################################

    
    
    

    

    cat('--------------- Organizing the results...')
    for(i in c(1:Boot_iters))
    {
      r <- result[[i]]

      SP <- rbind(SP,t(r[[1]]))
      SE <- rbind(SE,r[[2]])
    }
    cat('Done\n')
    #cat(paste0("N = ",as.character(N)," : Bootstrapped NMF done\n"))
    
    
    
    cat('--------------- Clustering the results\n')
    
    ## Step 5 (Clustering)
    incProgress(increment_progress_by, detail = paste0("N = ",as.character(N)," : Clustering"))
    
    
    ##############################################################################################################################
    ##############################################################################################################################
    ang_dis <- function(a,b)
    {
      if(((a%*%b)/(sqrt((a%*%a)*(b%*%b)))) > 1){return(0)}
      return((2*acos((a%*%b)/(sqrt((a%*%a)*(b%*%b)))))/pi)
    } 
    dist_from_cen <- function(i,j) return(ang_dis(SP[i,],centroids[j,]))       # return distance between object  and centroid j
    min_dist  <- function(i) return(which.min(mapply(dist_from_cen,i,c(1:N)))) # return the nearest centroid to the object i
    
    while(TRUE)
    {
      centroids <- matrix(0,1,K)        # intialize centroids using K-means++ algorithm
      cen_indxs <- ceiling(runif(1)*(dim(SP)[1]))    # choose the first centroid randomely
      rem_indxs <- setdiff(c(1:(dim(SP)[1])),cen_indxs)
      centroids[1,] <- SP[cen_indxs[1],]
      min_distance_from_centroids <- function(i) # return min distance of object i from currently selected centroids
      {
        dists <- sapply(c(1:(dim(centroids)[1])),function(j) ang_dis(SP[i,],centroids[j,]))
        return(min(dists))
      }
      if(N > 1)
      {
        for(i in 1:(N-1)) # repeat for N-1 times to complete the set of N initial centroids
        {
          D_probs <- sapply(rem_indxs,min_distance_from_centroids)
          new_indx <- sample(rem_indxs,1,prob=D_probs)
          centroids <- rbind(centroids,SP[new_indx,])
          cen_indxs <- c(cen_indxs,new_indx)
          rem_indxs <- setdiff(rem_indxs,new_indx)
        }
      }
      
      clst0 <- rep(0,dim(SP)[1])     #SP[i,] <-> clst0[i]  (clusters indices: 1,2,3,...,N)
      clst1 <- rep(0,dim(SP)[1])     #SP[i,] <-> clst1[i]  (clusters indices: 1,2,3,...,N)
      
      num_iter <- 0
      while(TRUE) 
      {
        
        if(number_of_cpu_cores > 1) {
          cat('-------------------------- Updating cluster assignments...')
          clst1 <- sapply(c(1:(dim(SP)[1])),min_dist)
          cat('Done\n')
        } else {
          cat('-------------------------- Updating cluster assignments...')
          clst1 <- sapply(c(1:(dim(SP)[1])),min_dist)
          cat('Done\n')
        }

        for(i in c(1:N))
        {
          members_indices <- which(clst1 == i)
          if(length(members_indices) != 0)
          {
            centroids[i,] <- colMeans(matrix(SP[members_indices,],length(members_indices))) # calculate centroids
          }
        }
        
        # cost function: sum-of-squares
        cost <- sum(sapply(c(1:(dim(SP)[1])),function(i) ang_dis(SP[i,],centroids[clst1[i],])))
        if(identical(clst1,clst0)){break}
        num_iter <- num_iter + 1
        clst0 <- clst1
      }
      
      if(dim(table(clst1)) == N) {break} # if there is any single centroid, run the clustering procedure again!
    }
    ##############################################################################################################################
    ##############################################################################################################################
    
    
    history_of_centroids[[epoch_boot]] <- centroids
    
    if(epoch_boot %% 2 == 0) {
      old_epoch <- epoch_boot / 2
      old_centroids <- history_of_centroids[[old_epoch]]
      candidates <- list()
      for(j in 1:N) {
        distances <- sapply(c(1:N),function(i){ang_dis(centroids[j,],old_centroids[i,])})
        candidates[[j]] <- which(distances < Boot_conv)
      }
      assignment <<- rep(0,N)
      find_assignment <- function(current_index) {
        current_candidates <- candidates[[current_index]]
        possibles <- current_candidates[which(assignment[current_candidates] == 0)]
        if(length(possibles) == 0) {
          return(FALSE) 
        } 
        if(current_index == N) { 
          assignment[possibles] <<- current_index
          return(TRUE) 
        } 
        for(p in possibles) {
          assignment[p] <<- current_index
          if( find_assignment(current_index+1) == T) { return(TRUE) }
          assignment[p] <<- 0
        }
        return(FALSE)
      }
      can_find <- find_assignment(1)
      if(can_find == TRUE) {
        cat(paste0('--------------- Angular distance between centroids of epoch ',as.character(old_epoch),
                               ' and ',as.character(epoch_boot),' is less than ',as.character(Boot_conv),'\n'))
        break()
      } else {
        cat(paste0('--------------- Angular distance between centroids of epoch ',as.character(old_epoch),
                   ' and ',as.character(epoch_boot),' is more than ',as.character(Boot_conv),'\n'))
      }
    }
    
  }
  
  
  
  
  
  cat(paste0('----- ','Evaluating the clusters and saving the final results\n'))
  
  ## Step 6 (Evaluate)
  incProgress(increment_progress_by, detail = paste0("N = ",as.character(N)," : Evaluation of clustering result..."))
  
  
  
  avg_dist_from_clst <- function(i,j)   # average distance between object i and all members of cluster j
  {
    members <- which(clst1==j)
    return((sum(sapply(members,function(t) ang_dis(SP[i,],SP[t,]))))/length(members))
  }
  silh_width <- function(i)  # silhouette width of cluster i
  {
    a <- avg_dist_from_clst(i,clst1[i])
    out_avgs <- sapply(setdiff(c(1:N),clst1[i]),function(t) avg_dist_from_clst(i,t))
    if(N == 1)
    {
      b <- 1
    } else {
      b <- min(out_avgs)  
    }
    return((b-a)/(max(a,b)))
  }
  avg_silh_widths <- mean(sapply(c(1:(dim(SP)[1])),silh_width))  # average silhouette width of the whole clustering
  
  exposures <- matrix(0,N,G)          # contruct matrix E according to clustering of P
  for(i in c(1:N)) {exposures[i,] <- colMeans(matrix(SE[which(clst1 == i),],length(which(clst1 == i))))}
  
  M.re <- t(centroids) %*% exposures # reconstruct M
  re.E <- norm((M - M.re),type = 'F')   # Reconstruction error
  
  # Normalize signatures by making each column of matrix P sum up to 1
  # and each row of the exposure matrix should be modified so that the value of 'P x E' remain unchanged
  for(indx in 1:N)
  {
    mul <- sum(centroids[indx,])
    centroids[indx,] <- centroids[indx,] / mul
    exposures[indx,] <- exposures[indx,] * mul
  }
  
  write(c(N,avg_silh_widths,re.E),file=paste0(destination_folder,"Evaluation.txt"),append=(N!=1))
  write.table(centroids,paste0(destination_folder,"P-n-",as.character(N),".txt"))
  if(length(zero_samples) != 0)
  {
    exposures_complete <- matrix(0,N,(G+length(zero_samples)))         # At the end, we include zero samples by inserting sxtra zero columns into matrix E.
    exposures_complete[,-zero_samples] <- exposures
    exposures <- exposures_complete
  }
  write.table(exposures,paste0(destination_folder,"E-n-",as.character(N),".txt"))
  
  cat('\n')
}














cat('\nAll calculations finished!\n')




########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

cat('\nPlotting the evaluation diagram...')

plot_eval_diagram <- function()
{
  e <- read.table(paste0(destination_folder,"Evaluation.txt"))
  n <- e[,1]
  repro <- e[,2]
  frobe <- e[,3]
  frobe <- frobe/max(frobe)
  
  ## add extra space to right margin of plot within frame
  par(mar=c(5, 5, 4, 6) + 0.5)
  
  ymin <- 0.1
  
  ## Plot the second plot and put axis scale on right
  plot(n, repro, pch=20, axes=FALSE, ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25), xlab="", ylab="",
       type="p",col="red", main="Evaluation for N",frame.plot = FALSE)
  
  grid(lwd = 2)
  
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  plot(n, frobe, pch=20,  xlab="", ylab="", ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25),
       axes=FALSE, type="p", col="blue",frame.plot = FALSE)
  ## a little farther out (line=4) to make room for labels
  mtext("Relative Frobenius Reconstruction Error",side=4,col="blue",line=2.5)
  axis(4,lwd = 2, ylim=range(frobe), col="blue",col.axis="blue",las=1)
  
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  ## Plot first set of data and draw its axis
  plot(n, repro, pch=20, axes=FALSE, ylim=c(0,1),xlim = c(n[1],n[length(n)]+0.25), xlab="", ylab="",
       type="p",col="red", main="Evaluation for N",frame.plot = FALSE)
  lines(n, repro,col='red')
  axis(2, lwd = 2,ylim=range(repro),col="red",col.axis="red",las=1)  ## las=1 makes horizontal labels
  mtext("Signatures Reproducibility",side=2,col ="red",line=3.75)
  
  ## Draw the time axis
  axis(1,(n[1]-1):(n[length(n)]+1))
  mtext("Number of mutational signatures",side=1,col="black",line=2.5)
}

pdf(paste0(destination_folder,"Evaluation_diagram.pdf"))
plot_eval_diagram()
dev.off()


e <- read.table(paste0(destination_folder,"Evaluation.txt"))
repro <- e$V2
Drop_in_repro <- sapply(c(1:(length(repro)-1)),function(i){return(repro[i]-repro[i+1])})
N_opt <- order(-Drop_in_repro)[1]
write.table(N_opt,file=paste0(destination_folder,"N_opt.txt"),row.names=F,col.names=F)

cat('Done\n')




