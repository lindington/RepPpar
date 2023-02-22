library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(RColorBrewer)

# Annotation file is in plink cluster format
# Read input file
setwd("C:/Users/Jag/Documents/LMU/work/01.capture/02.pca")
covar <- read.table('219_out_maxdepth.covar', stringsAsFactors = FALSE);
# Read annot file
annot <- read.table('plink_bait219.clst',sep="",header=TRUE);
#note that plink cluster files are usually tab - separated instead
#View(annot)

# Parse components to analyze
comp <- as.numeric(strsplit('1-2',"-",fixed=TRUE)[[1]])
#View(comp[1])
# Eigenvalues
eig <- eigen(covar,symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val,digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V","PC",colnames(PC))
PC$Pop <- factor(annot$CLUSTER)
Pops<-annot$CLUSTER[!duplicated(annot$CLUSTER)]


title <- paste("PC",comp[1],"(",signif(eig$val[comp[1]],digits=3)*100,"%)","/PC",comp[2],"(",signif(eig$val[comp[2]],digits=3)*100,"%)",sep="",collapse="")
#View(title)

plottit<-data.frame(PC$PC1,PC$PC2,PC$Pop,factor(annot$COLOUR))
colnames(plottit)<-c("PC1","PC2","Pop","colrs")

x_axis = paste("PC",comp[1],sep ="")
y_axis = paste("PC",comp[2],sep ="")

plottit$Pop<-factor(plottit$Pop, levels=Pops)

plottit$hz<-c(rep("Outgroup",10),rep("Portallet",120),rep("Pays Basco",89))
plottit$ancestry<-c(rep("parallelus",36),rep("hybrid",59),rep("erythropus",35),rep("parallelus",20),rep("hybrid",20),rep("parallelus",29),rep("hybrid",20))
colperpop<-c("black"="black","#A6CEED"="#94b2d6","#FAD9BE"="#f5bda4","#ED797C"="#e47c7c")
colperanc<-c("parallelus"="#94b2d6","hybrid"="#f5bda4","erythropus"="#e47c7c")

#plottit$colrs<-colperpop

ggplot(data = plottit, aes(x=PC1,y=PC2,fill=ancestry,shape=hz)) +
  geom_point(size=3,
             stroke=0.4)+
  scale_fill_manual(name = "Ancestry",
                    labels = c("erythropus", "hybrid", "parallelus"), 
                    breaks = waiver(), 
                    values = c("#e47c7c","#f5bda4","#94b2d6"))+
  scale_shape_manual(name = "Transect",
                     labels = c("Outgroup","Pays Basco", "Portalet"),
                     values = c(22,21,24))+
  scale_y_continuous(minor_breaks=NULL, limits=c(-0.16,0.08)) +
  scale_x_continuous(minor_breaks=NULL, limits=c(-0.09,0.11)) + 
  coord_fixed()  +
  theme_bw() +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle(title)

#                     values = c(15,16,17))


# ggplot() + 
#   geom_point(data=plottit, aes(x=PC1,y=PC2,fill=colrs),
#              colour="white",
#              pch = 21,  
#              size = 3, 
#              stroke = 0.4, 
#              show.legend = TRUE) +
#   scale_y_continuous(minor_breaks=NULL, limits=c(-0.16,0.1)) +
#   scale_x_continuous(minor_breaks=NULL, limits=c(-0.11,0.11)) + 
#   coord_fixed()  +
#   ggtitle(title) + 
#   scale_fill_manual(values=colperpop,name="NGSadmix K=3 assignment",breaks=c("black","#A6CEED","#FAD9BE","#ED797C"),labels=c("out","parallelus","hybrid", "erythropus"))
# 

ggsave('bait_PCA12.pdf')

unlink("Rplots.pdf",force=TRUE)

###################################
#### Start components 3 and 4: ####
###################################
comp <- as.numeric(strsplit('3-4',"-",fixed=TRUE)[[1]])
View(comp[1])

eig <- eigen(covar,symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val,digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V","PC",colnames(PC))
PC$Pop <- factor(annot$CLUSTER)
PC$Col <- factor(annot$COLOUR)
#View(PC)


title <- paste("PC",comp[1],"(",signif(eig$val[comp[1]],digits=3)*100,"%)","/PC",comp[2],"(",signif(eig$val[comp[2]],digits=3)*100,"%)",sep="",collapse="")
View(title)
x_axis = paste("PC",comp[1],sep ="")
y_axis = paste("PC",comp[2],sep ="")


ggplot() + 
  geom_point(aes_string(x=x_axis,y=y_axis, fill=PC$Pop),
             data=PC,
             colour="white",
             pch = 21,  
             size = 3, 
             stroke = 0.4, 
             show.legend = TRUE) +
  #scale_y_continuous(minor_breaks=NULL, limits=c(-0.2,0.2)) +
  #scale_x_continuous(minor_breaks=NULL, limits=c(-0.2,0.2)) + 
  #coord_fixed()  +
  ggtitle(title) + 
  scale_fill_manual(values=colperpop,name="Populations",breaks=c("GAB","HER","SOQ","TOU","ARA","POR","MUL","FOR","PAZ","LAN","ESC"),labels=c("Gabas","Hermine","Soques","Tourmont","Araille","Portalet","Mulas","Formigal","Pazino","Lanuza","Escarilla"))+
  
  ggsave('bait_PCA34.pdf')

unlink("Rplots.pdf",force=TRUE)
