pops = c("OUT","GAB","HER","SOQ","TOU","ARA","POR","MUL","TRO","FOR","PAZ","LAN","ESC","STP","AIN","URD","OTS","ARI","LEK","MAR","ALM","VEN")
#getwd()
setwd("/dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/10.pi_theta_taj")
Npop <- data.frame(10,5,10,11,10,10,13,13,13,10,10,10,5,10,10,10,10,10,9,10,10,10)
colnames(Npop) <- pops

for (i in 1:22){
    pop <- pops[i]
    print(pop)
    N <- as.numeric(Npop[i]*2) #number of haploid sequences
    theta = read.table(paste0("01.output/logThetas_",pop,".chr1"), header=T, comment.char="")
    genes = read.table("../00.input/bait.sites") #use gene position
    genes$V1 <- genes$V2
    genes$V2 <- genes$V3
    genes$V3 = NA
    genes$V4 = NA
#View(genes)
    
    startpos = genes$V1
    endpos = genes$V2
    genes$V3 <- as.vector(startpos)
    genes$V4 <- as.vector(endpos)


#coef needed for calculation of D
    a1 = sum(1/1:(N-1));  a2 = sum((1/1:(N-1))^2)
    b1 = (N+1)/(3*(N-1)); b2 = 2*(N^2 + N + 3)/(9*N*(N-1))
    c1 = b1 - 1/a1;       c2 = b2 - (N+2)/(a1*N) + a2/a1^2
    e1 = c1/a1;           e2 = c2/(a1^2+a2)

#make cols for D calcs
    genes$thetaW = NA
    genes$pi = NA
    genes$S = NA
    genes$TajimaD = NA
     
#loop over genes and calculate thetas, number variant sites, D
    curr_t<-1
    n_t<-nrow(theta)
    start_pos<-c()
    end_pos<-c()
    for(i in 1:nrow(genes)){
        not_discovered<-T
        for(cp in curr_t:n_t){
            if(not_discovered){
                if(theta$Pos[cp]>=genes$V3[i]){
                    start_pos<-c(start_pos,cp)
                    curr_t<-cp
                    not_discovered<-F
                }
            }
            if(theta$Pos[cp]>genes$V4[i]){
                end_pos<-c(end_pos,cp-1)
                break
            }
        }
        print(i)
    }
    
    if(length(start_pos)!=length(end_pos)){
        end_pos<-c(end_pos,nrow(theta))
    }
    
    sites = sapply(1:length(start_pos),function(i) start_pos[i]:end_pos[i])
    
    genes$thetaW = sapply(1:nrow(genes),function(i) sum(exp(theta$Watterson[sites[[i]]])))
    genes$pi = sapply(1:nrow(genes),function(i) sum(exp(theta$Pairwise[sites[[i]]])))
    
    genes$S = a1*genes$thetaW
    genes$TajimaD = (genes$pi - genes$thetaW) / sqrt(e1*genes$S + e2*genes$S^2)
    
    
#    View(genes)
#    View(data.frame(start_pos,end_pos))
    write.table(genes, file = paste0('01.output/per_gene_thetas_',pop),row.names=FALSE, quote = FALSE, sep = "\t")
    
}
