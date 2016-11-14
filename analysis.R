#===== ANALYSIS ======#

## Simulate an outbreak and reconstruct it using T, G and C ##
w.dens <- c(1,2,3,1,0.5,0.2)
outbreak <- outbreaker::simOutbreak(2,w.dens,n.hosts=30)
CTD <- clust.simCTD(outbreak,1,5)


min.eps <- 2*nrow(CTD)/((outbreak$n-2)*(outbreak$n+1)) 
max.eps <- nrow(CTD)/(outbreak$n-1)

prior.eps <- c(max(c(0,min.eps)),min(1,max.eps))

data=list(dna=NULL,dates=outbreak$onset,w.dens=w.dens,CTD=CTD)
config=list(n.iter=2e4, sample.every=200,find.import=FALSE)

result <- outbreaker(data=data, config=config)
analysis <- result.analysis(result,outbreak)


seq <- 3:30

lambda <- 1
mod1 <- function(n,lambda) 1/(1 + lambda*(1/2 - 1/n)) 
mod2 <- function(n,lambda) 1/(1+(lambda/2)*(n-2)/(n-1))
informative <- function(n,lambda) 1 - lambda/(n-2)
asymp <- 1/(1+lambda/2)
df <- data.frame(x=seq,mod1=mod1(seq,lambda),mod2=mod2(seq,lambda),inform=informative(seq,1),asymp=asymp)
df <- melt(df,id="x")

p <- ggplot(df,aes(x,value,colour=variable)) + geom_line()
p


## Will simulate outbreaks and CTD, and then plot prop(CTD==TP) against outbreak size, and plot theoretical value
plot.const.prop <- function(lambda,runs=100){
  
  if(lambda>40) break("lambda is too large")
  
  temp <- function(lambda,n){
    repeat({
      outbreak <- outbreaker::simOutbreak(2,w.dens,n.hosts=n)
      if(outbreak$n>(lambda+2)) break
    })
    CTD<-clust.simCTD(outbreak,1,lambda,full=TRUE)
    prop <- mean(sum(CTD$contact)/sum(CTD$accept))
    return(c(outbreak$n,prop))
  }
  
  out <- sapply(sample(40:60,runs,replace=TRUE),temp,lambda=lambda)
  plot(out[1,],out[2,])
  abline(h=1/(1+(lambda/2)))
  
}