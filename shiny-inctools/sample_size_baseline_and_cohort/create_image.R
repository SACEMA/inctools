create_plot <- function(samples, prevalence, reagents, type, min_pools, min_indiv){
  library(ggplot2)
#   prevalence=5
#   samples=100
#   reagents=7
#   type="Dried blood spot"
#   min_pools=20
#   min_indiv=10
  
  prevvec=seq(0,10,0.5)
  
  if(type=="Dried blood spot") possible_ps <- 2:5 else possible_ps <- 2:10
  
  opsom <- matrix(0,1,8)
  
  for(z in 1:length(prevvec)){
    
    vir_y <- matrix(0,length(possible_ps),8)
    
    for(y in 1:length(possible_ps)){
      prev <- prevvec[z]/100
      ps <- possible_ps[y]
      
      mat <- matrix(0, nrow=365, ncol=13)
      
      for(i in 1:365){
        if(i==1){
          p <- floor(samples/ps) #pools required to run all samples (but we do no put less samples into pool than specified!)
          x <- 0 #xtras due to positive pools from day before
          #change assumption: we will still run the pools to a minimum of 4 pools.
          #r <- floor( p/22)+ifelse(p/22-floor(p/22)<10/22, 0 , 1 ) #runs performed on day (when pooling)
          #l <-ifelse(p/22-floor(p/22)>=10/22, 0 ,samples-r*22*ps) #left over samples (when pooling)
          r <- floor(p/22)+ifelse(((p/22-floor(p/22))*22*ps)<min_pools, 0 , 1 )
          p_done <- floor(p/22)*22 + ifelse(((p/22-floor(p/22))*22)*ps<min_pools, 0 , ((p/22-floor(p/22))*22) ) 
          l <- samples - p_done*ps
          
          positives <- rbinom(1,samples,prev)
          positions <- sort(sample(samples,positives))
          pos_l <- positives - sum(positions<p_done*ps)
          positions <- positions[positions<p_done*ps]
          pos_pools <- length(unique(ceiling(positions/ps)))
          
          run_indiv <- floor(samples/22)+ifelse((samples/22-floor(samples/22))*22<min_indiv, 0 , 1 ) #runs performed on day (when not pooling)
          l_indiv <- ifelse((samples/22-floor(samples/22))*22>=min_indiv, 0 ,samples-run_indiv*22) #left over samples (when not pooling)
          
          
          cons_indiv <- ifelse(reagents==0, samples, samples*reagents)
        }
        else{
          samples_i <- samples+mat[i-1,6]
          positives <- rbinom(1,samples,prev) + mat[i-1,8]
          p <- floor(samples_i/ps) 
          #max number of samples you have to run again because of positive pools (assuming only 1 positive per pool)
          #x = ceiling(mat[i-1,6]*ps * prev)*ps
          #expected number of samples you have to run again because of positive pools (from yesterday)
          x = mat[i-1,9] *ps      
          
          r <- floor((p+x)/22)+ifelse((((p+x)/22-floor((p+x)/22))*22*ps)<min_pools, 0 , 1 )
          p_done <- floor((p+x)/22)*22 + ifelse((((p+x)/22-floor((p+x)/22))*22)*ps<min_pools, 0 , (((p+x)/22-floor((p+x)/22))*22) ) -x
          l <- round(samples_i - p_done*ps, 0)
          
          positions <- sort(sample(samples_i,positives))
          pos_l <- positives - sum(positions<p_done*ps)
          positions <- positions[positions<p_done*ps]
          pos_pools <- length(unique(ceiling(positions/ps)))
          
          samples_ind <- samples+mat[i-1,11]
          run_indiv <- floor( samples_ind/22)+ifelse((samples_ind/22-floor(samples_ind/22))*22<min_indiv, 0 , 1 )
          l_indiv <- ifelse((samples_ind/22-floor(samples_ind/22))*22>=min_indiv, 0 ,samples_ind-run_indiv*22)
          
          
          cons_indiv <- ifelse(reagents==0, (samples_ind-l_indiv), (samples_ind-l_indiv)*reagents)
        }
        
        cons_pool <- ifelse(reagents==0, (p_done+x), (p_done+x)*reagents) 
        cons_indiv <- cons_indiv
        
        mat[i,1] <- samples; mat[i,2] <- positives; mat[i,3] <- p ; mat[i,4] <- x ;  mat[i,5] <- r ; mat[i,6] <- l ; mat[i,7]<- p_done ; 
        mat[i,8] <- pos_l; mat[i,9] <- pos_pools;
        mat[i,10] <- run_indiv ; mat[i,11] <- l_indiv 
        mat[i,12] <- cons_pool ;  mat[i,13] <- cons_indiv
        
      }
      
      colnames(mat) <- c("samples","positives" ,"num_pools","xtras","runs","samples_left", "pools_done","positives_left","positive pools" ,"runs_indiv",
                         "samples_left_indiv", "cons_pool","cons_indiv")
      mat <-  as.data.frame(mat)
      
      
      runs_pool <- mat$runs
      runs_indiv <- mat$runs_indiv
      runs_saved <- round((1- sum(runs_pool)/sum(runs_indiv))*100,2) 
      
      cost_pool <- sum(mat$cons_pool)
      cost_indiv <- sum(mat$cons_indiv ) 
      perc_saved <-   round((1-((1/ps)+(1-dbinom(0,ps,prev))))*100,2) #round((1- cost_pool/cost_indiv)*100,2) 
      
      vir_y[y,] <- c(prev*100, ps,median(runs_pool),median(runs_indiv),runs_saved, cost_pool/365, cost_indiv/365, perc_saved)
      
    }
    opsom <- rbind(opsom, vir_y)
  }
  colnames(opsom) <- c("prevalence","poolsize","runs_pool","runs_indiv","runs_saved",
                       "cost_pool", "cost_indiv","cost_saved" )
  opsom <- as.data.frame(opsom[-1,])
  opsom$prevalence=rep(prevvec, each=length(possible_ps))
  opsom$cost <- 100-opsom$cost_saved
  

  long <- opsom[,c(1,2,9)]
  library(tidyr)
  short <- spread(long,prevalence,cost)
  short <- t(as.matrix(short))
  short <- short[-1,]
  min_cost <- opsom$cost[which(opsom$prevalence==prevalence & opsom$cost==min(opsom$cost[opsom$prevalence==prevalence]))]
  ps <- max(opsom$poolsize[which(opsom$prevalence==prevalence & opsom$cost==min(opsom$cost[opsom$prevalence==prevalence]))])
  
  amount<- round(opsom$cost_indiv[opsom$prevalence==prevalence & opsom$poolsize==ps]-
                   opsom$cost_pool[opsom$prevalence==prevalence & opsom$poolsize==ps],0)
  
  percentage<- round(opsom$cost_saved[opsom$prevalence==prevalence & opsom$poolsize==ps],0)
  
  runs_text1 <- opsom$runs_indiv[opsom$prevalence==prevalence & opsom$poolsize==ps]
  runs_text2 <- opsom$runs_pool[opsom$prevalence==prevalence & opsom$poolsize==ps]
  
  if(runs_text1>runs_text2){  
    if(reagents>0) {
      cost_text <- paste("Based on your input values, at the optimal pool size of ",ps, ", you will save ", percentage, "% of costs and daily batched
                                    runs will decrease from ", runs_text1, " to ", runs_text2, ". Your estimated daily cost saving is $", amount,".", sep="")
      
    } else { 
      cost_text <- paste("Based on your input values, at the optimal pool size of ",ps, ", you will save ", percentage, "% of costs and daily batched
                                    runs will decrease from ", runs_text1, " to ", runs_text2, ".", sep="")
    }
  }  else {
    if(reagents>0) {
      cost_text <- paste("Based on your input values, at the optimal pool size of ",ps, ", you will save ", percentage, "% of costs. Your estimated daily cost saving is $", amount,".", sep="")
    }   else {
      cost_text <- paste("Based on your input values, at the optimal pool size of ",ps, ", you will save ", percentage, "% of costs.", sep="")
    }
  }

  
  
  
  
  return(list(cost_text=cost_text,prevvec=prevvec, possible_ps=possible_ps, ps=ps, min_cost=min_cost, short=short, opsom=opsom))
}