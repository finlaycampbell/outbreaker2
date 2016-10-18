#===== ANALYSIS ======#

## Simulate an outbreak and reconstruct it using T, G and C ##

w.dens <- c(1,2,3,1,0.5,0.2)
outbreak <- outbreaker::simOutbreak(2,w.dens,n.hosts=30)
CTD <- simCTD(outbreak,0.3,0.8)

result <- outbreaker(data=list(dna=NULL,dates=outbreak$onset,w.dens=w.dens,CTD=NULL),
                                  config=list(n.iter=2e4, sample.every=200,find.import=FALSE))

result2 <- outbreaker(data=list(dna=NULL,dates=outbreak$onset,w.dens=w.dens,CTD=CTD),
                     config=list(n.iter=2e4, sample.every=200,find.import=FALSE))

analysis <- result.analysis(result,outbreak)
analysis2 <- result.analysis(result2,outbreak)

#===== FUNCTIONS ======#

## Return the likelihood of two individuals being a transmission pair depending ##
## on their contact status, given eps and xi                                    ##
contact.ll <- function(eps,xi,cij){
  
  if(!inherits(cij,"logical")) stop("cij is not a logical")
  
  #If the coverage is 0, contacts cannot inform us of the likelihood of a given transmision pair
  if(cij) eps/(eps+eps*xi)
  
  #If the coverage is 1, the probability of being a tranmission without a contact is zero
  else if(!cij) (1-eps)/(2-eps-eps*xi)
  
}

## Plot the support provided for ancestries with reported contacts, across eps and xi ##
plot.support <- function(){
  
  eps.seq <- seq(0.1,0.9,length=20)
  xi.seq <- seq(0,2,length=20)
  out <- expand.grid(eps=eps.seq,xi=xi.seq)
  
  out$weight <- apply(out,1,function(x) contact.ll(x[1],x[2],TRUE)/contact.ll(x[1],x[2],FALSE)-1)
  
  p <- ggplot(out,aes(out$eps,out$xi,fill=weight)) + geom_tile() + 
       scale_fill_gradient2(low="firebrick",mid="white","high"="darkgreen")
  
  return(p)

}

## A function to analyse the outbreaker output for accuracy, entropy and confidence 
result.analysis <- function(result,true.outbreak,plot=TRUE,print=FALSE){
  id <- seq_len(true.outbreak$n)
  adder <- which(names(result)=="alpha.1")-1
  samples <- length(result$step)
  
  #Determine the modal transmission network
  network <- data.frame(from=do.call(rbind,lapply(id, function(i) ({
    modal.ances <- as.integer(names(which.max(table(result[[i+adder]]))))
    if(length(modal.ances)==0) return(NA) else return(modal.ances)
  }))), to=id)
  
  import <- which(is.na(network$from))
  
  #Define the indices to call the times of infection
  t.inf.index <- which(names(result)=="t.inf.1")-1+id
  
  #Determine the median posterior time of onset for plotting
  onset <- unlist(lapply(result[t.inf.index],median))
  
  #Scale onset to begin at 0
  onset <- onset - min(onset)
  
  #Determine confidence in our results
  
  transmission.id <- id[!sapply(result[id+adder],function(i) any(is.na(i)))]
  
  confidence <- round(mean(sapply(transmission.id,function(i) mean(result[[i+adder]]==true.outbreak$ances[i])),na.rm=TRUE),2)
  
  entropy <- round(mean(sapply(result[transmission.id+adder],function(i) {fk <- table(i)/sum(table(i))
  -sum(log(fk)*fk)})),2)
  
  #Determine the proportion of correctly inferred ancestries
  num.correct <-  sum(true.outbreak$ances==network$from,na.rm=TRUE)
  num.correct <- num.correct + sum(is.na(true.outbreak$ances[is.na(network$from)]))
  accuracy <- round(num.correct/nrow(network),2)
  
  out <- list(analysis=data.frame(accuracy=accuracy,confidence=confidence,entropy=entropy))
  
  return(out)
}