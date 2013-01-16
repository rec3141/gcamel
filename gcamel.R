rm(list=ls())
# parameters needed:
#  1) GC content of genome at each position (.ffn file)
#  2) GC content of xene at each position (.fasta file)
#  3) transition/transversion ratio of genome (assume 2:1)
#  4) substitution rate of genome (assume 0.0245% and 0.455% per Myr for nonsynonymous and synonymous)
#  5) linear (?) relationships between {GC1,GC2} and {GC3} from Muta & Osawa
# 
# wanted:
#  1) GC content of xene at introduction
#  2) time since introduction of xene

# method to get 1) and 2)
#  a) calculate delGC for 1Myr time interval
#     delGC = S * (tstv+0.5)/(tstv+1) * (gc_genome - gc_xene)
#  b) add delGC to gc_xene
#  c) calculate correlation between GC1,GC2~GC3 and Muta&Osawa
#  d) repeat a-c for 1000Myr
#  e) find minimum correlation for all c) -- this point is initial GC and time since transfer


fit.sigmoid.gc1 <- function(x_in) {
#zunzun.com
#finds gc3 given gc1

	a <- 1.4180794193642083E+00
	b <- 1.2403719287318339E+02
	c <- 7.5494898937434289E+00

	y = a / (1.0 + b*exp(-1.0 * c * x_in))
}

fit.sigmoid.gc2 <- function(x_in) {
#zunzun.com
#finds gc3 given gc2

	a = 1.0734349192542731E+00;
	b = 2.7879639168098345E+02;
	c = 1.3906208449024092E+01;

	y = a / (1.0 + b*exp(-1.0 * c * x_in));
}

reps <- 1000 #bootstrap reps
mya <- 1000 #millions of years to ameliorate
tstv <- 2
Srate <- c(0.123,0.045,0.668)/100/2 #weighted substitution rates

file.genome <- 'genome.ffn'
gcs.host <- as.numeric(unlist(strsplit(system("perl gc_by_pos.pl 0 1 genome.ffn",intern=TRUE),"\t")))
gcs.host.genes <- do.call(cbind, strsplit(system("perl gc_by_pos.pl 1 1 genome.ffn",intern=TRUE),"\t"))

genome.gcs <- read.table(file="genome_gcs.csv")
list.xenes <- list.files(pattern="soxBDAG.*fst")
tminmax <- NULL
oldgc <- NULL
# closest.genome <- NULL

#   par(mfrow = c(2, 1))
  xs <- 1:100/100
  gc1.fit <- fit.sigmoid.gc1(xs)
  gc2.fit <- fit.sigmoid.gc2(xs)

#   points(genome.gcs[,1],genome.gcs[,3],pch='.',col='grey')
#   points(genome.gcs[,2],genome.gcs[,3],pch='.',col='grey')
# 
#   points(gcs.host.genes[1,],gcs.host.genes[3,],pch='.',col='pink')
#   points(gcs.host.genes[2,],gcs.host.genes[3,],pch='.',col='purple')

for (myfile in list.xenes) {
  file.xene <- myfile
  print(file.xene)
  plot(xs,gc1.fit,xlim=c(0.2,0.8),ylim=c(0,1),xlab='GC1,GC2',ylab='GC3',type='l',col='pink',main=file.xene)
  lines(xs,gc2.fit,col='purple')

  gcs.xene.get <- as.numeric(unlist(strsplit(system(paste("perl gc_by_pos.pl 0 ",reps,file.xene,sep=" "),intern=TRUE),"\t")))
  gcs.xene <- matrix(gcs.xene.get,ncol=3,byrow=TRUE)

  points(gcs.host[1],gcs.host[3],pch=3,col='red')
  points(gcs.host[2],gcs.host[3],pch=3,col='red')

  points(gcs.xene[,1],gcs.xene[,3],pch='.',col='pink')
  points(gcs.xene[,2],gcs.xene[,3],pch='.',col='purple')

 for (rep in 1:reps) {
  print(rep);
  rev.xene.save <- NULL
  fwd.xene.save <- NULL

  delGC.fwd <-c(0,0,0)
  delGC.rev <-c(0,0,0)

  fwd.xene <- gcs.xene[rep,]
  rev.xene <- gcs.xene[rep,]

#   points(gcs.xene[rep,1],gcs.xene[rep,3],pch='o',col='pink')
#   points(gcs.xene[rep,2],gcs.xene[rep,3],pch='o',col='purple')

  dist.gc1 <- rep.int(0,mya)
  dist.gc2 <- rep.int(0,mya)
  dist.both <- rep.int(0,mya)

  for (i in 1:mya) {
    delGC.fwd <- Srate * (tstv+0.5)/(tstv+1) * (gcs.host - fwd.xene)
    delGC.rev <- Srate * (tstv+0.5)/(tstv+1) * (gcs.host - rev.xene)

    #forward
    fwd.xene <- fwd.xene + delGC.fwd;
    fwd.xene.save <- rbind(fwd.xene.save,fwd.xene)
    #reverse
    rev.xene <- rev.xene - delGC.rev;
    rev.xene.save <- rbind(rev.xene.save,rev.xene)

#     genome.dists <- rowSums((genome.gcs - rev.xene)^2)
#     close.genome <- rbind(close.genome,c(which.min(genome.dists),min(genome.dists)))

  }
  lines(fwd.xene.save[,1],fwd.xene.save[,3],col='red')
  lines(fwd.xene.save[,2],fwd.xene.save[,3],col='red')
  lines(rev.xene.save[,1],rev.xene.save[,3],col='pink')
  lines(rev.xene.save[,2],rev.xene.save[,3],col='purple')

  dist.gc1 <- abs(fit.sigmoid.gc1(rev.xene.save[,1]) - rev.xene.save[,3])
  dist.gc2 <- abs(fit.sigmoid.gc2(rev.xene.save[,2]) - rev.xene.save[,3])
  dist.both <- (fit.sigmoid.gc1(rev.xene.save[,1]) - rev.xene.save[,3])^2 + (fit.sigmoid.gc2(rev.xene.save[,2]) - rev.xene.save[,3])^2

  tminmax <- rbind(tminmax,c(file.xene,rep,which.min(dist.gc2),which.min(dist.gc1),which.min(dist.both)))
  oldgc <- rbind(oldgc,c(file.xene,rep,rev.xene.save[which.min(dist.gc2)],rev.xene.save[which.min(dist.gc1)],rev.xene.save[which.min(dist.both)]))
#   closest.genome <- rbind( closest.genome, close.genome[ which.min(close.genome[,2]) , ] )

 }

#   points(gcs.host[1],gcs.host[3],pch=3,col='red')
#   points(gcs.host[2],gcs.host[3],pch=3,col='red')

#   xs <- 1:reps
#   plot(dist.gc1~xs,col='violet',type='l',ylim=c(0,1))
#   lines(dist.gc2~xs,col='blue')
#   lines(dist.both~xs,col='green')

}

tempfiles <- tempfile("gcamel-",tmpdir=".")
write.csv(tminmax,file=paste(tempfiles[1],"-tminmax.csv",sep=""))
write.csv(oldgc,file=paste(tempfiles[1],"-oldgc.csv",sep=""))
# write.csv(closest.genome,file=paste(tempfiles[1],"-closest.csv",sep=""))
