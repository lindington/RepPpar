#Rscript -e 'write.table(cbind(seq(1,50),rep(1,50),c(rep("ARA",10),rep("ESC",5),rep("FOR",8),rep("GAB",5),rep("HER",10),rep("LAN",10),rep("MUL",7),rep("SOQ",1),rep("MUL",2),rep("PAZ",10),rep("POR",10),rep("SOQ",10),rep("TOU",10))),row.names=F,col.names=c("FID","IID","CLUSTER"),sep ="\t",file="plink_bait98.clst",quote=F)'

# seq & rep 1,50 or 1,N?

Rscript -e 'write.table(cbind(seq(1,219),rep(1,219),c(rep("OUT",10),rep("GAB",5),rep("HER",10),rep("SOQ",11),rep("TOU",10),rep("ARA",10),rep("POR",13),rep("MUL",13),rep("TRO",13),rep("FOR",10),rep("PAZ",10),rep("LAN",10),rep("ESC",5),rep("STP",10),rep("AIN",10),rep("URD",10),rep("OTS",10),rep("ARI",10),rep("LEK",9),rep("MAR",10),rep("ALM",10),rep("VEN",10))),row.names=F,col.names=c("FID","IID","CLUSTER"),sep ="\t",file="plink_bait219.clst",quote=F)'
