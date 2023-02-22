install.packages("data.table")
install.packages("ggplot2")
install.packages("plyr")
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GeneOverlap")
install.packages("overlapping")


library(data.table)
library(ggplot2)
library(plyr)
library(tidyverse)
library(GeneOverlap)
library(boot)
library(overlapping)

getwd()

## finding correlation between genes with high fst and genes with high dxy
#read in fst values for each pop
por_fst <- fread("genefst_GAB_ESC", header = FALSE)
bas_fst <- fread("genefst_STP_VEN", header = FALSE, sep = " ")
colnames(por_fst) <- c("V1","V2","V3","fst","siteswithdata")
#read in dxy values
por_dxy <- fread("03.dxy/03.output/per_gene_dxy_por.csv", header = TRUE)
por_dxy <- por_dxy %>%
  filter(por_dxy$siteswithdata > 1000)
bas_dxy <- fread("03.dxy/03.output/per_gene_dxy_bas.csv", header = TRUE)
bas_dxy <- bas_dxy %>%
  filter(bas_dxy$siteswithdata > 1000)

# only keep genes that exist in both sumstats for both hz:
por_fstdxy <- por_dxy$V1[(por_dxy$V1 %in% por_fst$V2)]
bas_fstdxy <- bas_dxy$V1[(bas_dxy$V1 %in% bas_fst$V2)]

hzs_fstdxy <- bas_dxy$V1[(bas_dxy$V1 %in% por_dxy$V1)]

por_dxy<-por_dxy %>%
  filter(V1 %in% por_fstdxy) %>%
  filter(V1 %in% hzs_fstdxy)
por_fst<-por_fst %>%
  filter(V2 %in% por_fstdxy) %>%
  filter(V2 %in% hzs_fstdxy)


bas_dxy<-bas_dxy %>%
  filter(V1 %in% bas_fstdxy) %>%
  filter(V1 %in% hzs_fstdxy)
bas_fst<-bas_fst %>%
  filter(V2 %in% bas_fstdxy) %>%
  filter(V2 %in% hzs_fstdxy)

# making pretty plots
por_genes <- cbind.data.frame(por_dxy$V1,por_dxy$dxy,por_fst$V5)
bas_genes <- cbind.data.frame(bas_dxy$V1,bas_dxy$dxy,bas_fst$V5)
colnames(por_genes) <- c("gene","dxy","fst")
colnames(bas_genes) <- c("gene","dxy","fst")

genes<-cbind.data.frame(por_genes$gene,por_genes$dxy,por_genes$fst,bas_genes$dxy,bas_genes$fst)
#order: genes por_dxy por_fst bas_dxy bas_fst
colnames(genes)<-c("gene","por","por_fst","bas","bas_fst")
#View(genes)


dxys<-stack(genes,select = c(por,bas),drop = TRUE)
colnames(dxys)<-c("dxy","hz")
#view(dxys)
fsts<-stack(genes,select = c(por_fst,bas_fst),drop = TRUE)
colnames(fsts)<-c("fst","hz")
#view(fsts)
loci<-rep(genes$gene,2)
#view(loci)

plotgenes<-cbind.data.frame(loci,dxys$dxy,fsts$fst,dxys$hz)
colnames(plotgenes)<-c("loci","dxy","fst","hz")

ggplot(data=plotgenes, aes(x=fst,y=dxy,color=hz))+
  geom_point(aes(x=fst,y=dxy,color=hz),size=1,shape=21)+
  geom_smooth(aes(fill=hz),method = lm,se=FALSE,size=1.2)+
  scale_color_manual(values = c("#D10000","#489CDA"),name="Hybrid Zone",label=c("porenees","bass"))+
  scale_fill_manual(values = c("#930D01","#26549C"),name="Hybrid Zone",label=c("porenees","bass"))+
  labs(title=expression("Per-gene d"["XY"]*"~F"["ST"]*" correlation"),x=expression("F"["ST"]),y=expression("d"["XY"]))+
  theme_classic()

#correlation in the hzs? linear regression dxy response, fst fixed

reg1 <- lm(por_genes$dxy~por_genes$fst,data=por_genes) 
summary(reg1)

  # Call:
  #   lm(formula = por_genes$dxy ~ por_genes$fst, data = por_genes)
  # 
  # Residuals:
  #   Min        1Q    Median        3Q       Max 
  # -0.013091 -0.003573 -0.001170  0.001922  0.032475 
  # 
  # Coefficients:
  #   Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    0.0148754  0.0001881   79.09   <2e-16 ***
  #   por_genes$fst -0.0150551  0.0004205  -35.80   <2e-16 ***
  #   ---
  #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  # 
  # Residual standard error: 0.005611 on 5581 degrees of freedom
  # Multiple R-squared:  0.1867,	Adjusted R-squared:  0.1866 
  # F-statistic:  1282 on 1 and 5581 DF,  p-value: < 2.2e-16

hist(residuals(reg1),
     col="darkgray")
plot(fitted(reg1),
     residuals(reg1))

reg2 <- lm(bas_genes$dxy~bas_genes$fst,data=bas_genes)
summary(reg2)

  # Call:
  #   lm(formula = bas_genes$dxy ~ bas_genes$fst, data = bas_genes)
  # 
  # Residuals:
  #   Min        1Q    Median        3Q       Max 
  # -0.012419 -0.004314 -0.001603  0.002197  0.040876 
  # 
  # Coefficients:
  #   Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)    0.0135781  0.0001763   77.03   <2e-16 ***
  #   bas_genes$fst -0.0133129  0.0005186  -25.67   <2e-16 ***
  #   ---
  #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  # 
  # Residual standard error: 0.006624 on 5581 degrees of freedom
  # Multiple R-squared:  0.1056,	Adjusted R-squared:  0.1054 
  # F-statistic:   659 on 1 and 5581 DF,  p-value: < 2.2e-16

hist(residuals(reg2),
     col="darkgray")
plot(fitted(reg2),
     residuals(reg2))
