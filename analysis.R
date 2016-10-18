#===== analysis ======#

w.dens <- c(1,2,3,1,0.5,0.2)
outbreak <- outbreaker::simOutbreak(2,w.dens)
CTD <- simCTD(outbreak,0.8,0.2)

result <- outbreaker2::outbreaker(data=list(dna=outbreak$dna,dates=outbreak$onset,w.dens=w.dens,CTD=CTD),
                                  config=list(n.iter=2e4, sample.every=200))