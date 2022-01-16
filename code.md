#Plotting code for forest plot（Figrue2A and Figure2B）
library(survival)
library(forestplot)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
#Univariate cox and multivariate cox regression analysis
#read data
tcga.expr <- read.table("tcga.expr.txt", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
tcga.surv <- read.table("tcga.surv.txt", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
#Keep samples with both expression data and survival data
intersect.tcga <- intersect(rownames(tcga.surv), colnames(tcga.expr))
tcga.expr <- tcga.expr[,intersect.tcga]
tcga.surv <- tcga.surv[intersect.tcga,]
#univariate cox
unicox <- data.frame()
for(i in 1:nrow(tcga.expr)){
  
  display.progress(index = i, totalN = nrow(tcga.expr))
  gene <- rownames(tcga.expr)[i]
  tmp <- data.frame(expr = as.numeric(tcga.expr[i,]),
                    futime = tcga.surv$OS.time,
                    fustat = tcga.surv$OS,
                    stringsAsFactors = F)
  cox <- coxph(Surv(futime, fustat) ~ expr, data = tmp)
  coxSummary <- summary(cox)
  unicox <- rbind.data.frame(unicox,
                             data.frame(gene = gene,
                                        HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                        z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                        pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                        lower = as.numeric(coxSummary$conf.int[,3][1]),
                                        upper = as.numeric(coxSummary$conf.int[,4][1]),
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}
#multivariate cox
mulcox.tcga <- summary(coxph(Surv(OS.time, OS) ~ ., data = tcga.surv))
mulcox.tcga <- data.frame(variable = rownames(mulcox.tcga$conf.int),
                          HR = mulcox.tcga$conf.int[,1],
                          lower.95CI = mulcox.tcga$conf.int[,3],
                          upper.95CI = mulcox.tcga$conf.int[,4],
                          p = mulcox.tcga$coefficients[,5],
                          stringsAsFactors = F)
#forestplot
tabletext[2,] <- c("Unicox",NA,NA,NA,NA) 
tabletext[9,] <- c("Multicox",NA,NA,NA,NA) 
pdf("forestplot of risk table.pdf", width = 8, height = 5)
forestplot(labeltext=tabletext,
           mean=c(NA,log2(as.numeric(hrtable$HR))),
           lower=c(NA,log2(as.numeric(hrtable$lower.95CI))), 
           upper=c(NA,log2(as.numeric(hrtable$upper.95CI))),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=0,
           lwd.zero=2,
           xticks = c(0,1,2,3,4),
           lwd.xaxis=2,
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2),
                           "9" = gpar(lwd=1, col="grey50", lty=2),
                           "12" = gpar(lwd=2, col="black")),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(.75,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())
pdf("forestplot of risk table.pdf", width = 8, height = 5)
forestplot(labeltext=tabletextICI,
           mean=c(NA,cox$OR),
           lower=c(NA,cox$OR.95L), 
           upper=c(NA,cox$OR.95H),
           graph.pos=6,
           graphwidth = unit(.25,"npc"),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),
           boxsize=0.4,
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=F,
           zero=1,
           lwd.zero=1,
           xticks = c(-1,0,1,2,3,4,5,6,7),
           lwd.xaxis=2,
           xlab=expression("log"[2]~"OR"),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),
                           "2" = gpar(lwd=1, col="grey50", lty=2)
                          ),
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           lineheight = unit(.75,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
invisible(dev.off())
#KM plot(Figure2C)
library("survival")
library("survminer")
fit <- survfit(Surv(time, status) ~ sex, data = survival)
summary(fit)
ggsurvplot(fit,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, 
          risk.table.col = "strata", 
          linetype = "strata", 
          surv.median.line = "hv", 
          ggtheme = theme_bw(), 
          palette = c("#E7B800", "#2E9FDF"))
#Sankey diagram（Figure2D and Figure2E）
#CAMOIP was used to make volcano plot(Figure2F)
#Plotting code for Figure3A
library(TCGAbiolinks)
library(maftools)
clinical <- GDCquery(project = "TCGA-STAD", 
                  data.category = "Clinical", 
                  file.type = "xml")

GDCdownload(clinical)

cliquery <- GDCprepare_clinic(clinical, clinical.info = "patient")
colnames(cliquery)[1] <- "Tumor_Sample_Barcode"
mut <- GDCquery_Maf(tumor = "STAD", pipelines = "mutect2")
maf <- read.maf(maf = mut, clinicalData = cliquery, isTCGA = T)
#plotting
col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site')
racecolors = RColorBrewer::brewer.pal(n = 4,name = 'Spectral')
names(racecolors) = c("ASIAN", "WHITE", "BLACK OR AFRICAN AMERICAN",  "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER")
gendercolors = c("SandyBrown","CadetBlue")
names(gendercolors) = c("FEMALE","MALE")
MUCcolors = c("SandyBrown","CadetBlue")
names(MUCcolors) = c("MT","WT")
ethnicitycolors=c("LightYellow1","LightYellow2")
names(ethnicitycolors)=c("HISPANIC OR LATINO","NOT HISPANIC OR LATINO")
gradecolors=c("LightSeaGreen","DarkOliveGreen1","SkyBlue","Goldenrod2")
names(gradecolors)=c("G1","G2","G3","GX")
annocolors = list(MUC16=MUCcolors,race_list = racecolors, 
                  gender = gendercolors, 
                  ethnicity = ethnicitycolors, 
                  neoplasm_histologic_grade = gradecolors
                  )
pdf("oncoplotTop20_col.pdf",width = 15,height = 10)
oncoplot(maf = maf,
         colors = col,
         annotationColor = annocolors,
         top = 20,
         clinicalFeatures = c("MUC16","race_list","gender","ethnicity","neoplasm_histologic_grade"),
         sortByAnnotation = TRUE,
         writeMatrix =T)
dev.off()
#plotting Figure3B
col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Intron')

typecolors = RColorBrewer::brewer.pal(n = 3,name = 'Spectral')
names(typecolors) = c("CIN", "GS", "EBV")
pathcolors = RColorBrewer::brewer.pal(n = 7,name = 'Spectral')
names(pathcolors) = c("Adenoneuroendocrine","M/D adeno ","Mixed (W/D and P/D)","P/D adeno","P/D adeno with NE feature","Signet ring cell ","W/D adeno")
rescolors=c("LightSeaGreen","DarkOliveGreen1","SkyBlue","Goldenrod2")
names(rescolors)=c("CR","PD","PR","SD")
MUCcolors=c("LightSeaGreen","DarkOliveGreen1")
names(MUCcolors)=c("WT","MT")
annocolors = list(TCGA = typecolors, 
                  Pathology = pathcolors, 
                  Best_of_response = rescolors,
                  MUC16=MUCcolors)
                  
pdf("oncoplotTop20_ICIcol.pdf",width = 15,height = 10)
oncoplot(maf = laml,
         colors = col,
         annotationColor = annocolors, 
         top = 20,
         clinicalFeatures = c("MUC16","TCGA","Pathology","Best_of_response"),
         sortByAnnotation = TRUE,
         writeMatrix =T)
dev.off()
#Plotting code for Figure3C and Figure3D
#Figure3C
output <- somaticInteractions(maf=STAD, top=50, pvalue=c(0.05, 0.01))
write.table(output$pairs, file="somaticInteractions.pairwise.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(output$gene_sets, file="somaticInteractions.comet.tsv", quote=FALSE, row.names=FALSE, sep="\t")
#Figure3D
output_ICI <- somaticInteractions(maf=maf, top=50, pvalue=c(0.05, 0.01))
write.table(output$pairs, file="somaticInteractions.pairwise.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(output$gene_sets, file="somaticInteractions.comet.tsv", quote=FALSE, row.names=FALSE, sep="\t")
#Plotting code for Figure3E
library(maftools)
library(trackViewer)
library(RColorBrewer)
df<-read.table("lolliplot.txt",as.is = T,sep = "\t",header = T)
newdf<-df[!is.na(df$HGVSp),]
gff = readRDS(file = system.file('extdata', 'protein_domains.RDs', package = 'maftools'))
refseqid <- strsplit(x = as.character(newdf$RefSeq), split = '.', fixed = TRUE)[[1]][1]
protein_inform <- gff[HGNC %in% "MUC16"][refseq.ID == refseqid,]
extractpos <- function(maf_aachange){
  prot.spl = strsplit(x = as.character(maf_aachange), split = '.', fixed = TRUE)
  prot.conv = sapply(sapply(prot.spl, function(x) x[length(x)]), '[', 1)
  pos = gsub(pattern = 'Ter.*', replacement = '',x = prot.conv)
  pos = gsub(pattern = '[[:alpha:]]', replacement = '', x = pos)
  pos = gsub(pattern = '\\*$', replacement = '', x = pos)
  pos = gsub(pattern = '^\\*', replacement = '', x = pos)
  pos = gsub(pattern = '\\*.*', replacement = '', x = pos)
  pos = as.numeric(sapply(X = strsplit(x = pos, split = '_', fixed = TRUE), FUN = function(x) x[1])) 
  aa = paste0(unlist(regmatches(maf_aachange, gregexpr("p[.].[0-9]+",maf_aachange))),"X")
  mutpos = data.frame(position = pos, mutation = maf_aachange, aa = aa, stringsAsFactors = F)
  return(mutpos[order(mutpos$pos),])
}
pos <- extractpos(newdf$HGVSp_Short)
head(pos)
nrpos <- pos[!duplicated(pos),]
rownames(nrpos) <- nrpos$mutation
nrpos$rate <- 1
nrpos[names(table(pos$mutation)),]$rate <-table(pos$mutation)
head(nrpos)
features <- GRanges("chr1", IRanges(start=protein_inform$Start,end=protein_inform$End,names=protein_inform$Description))
features$height <- 0.07 
features$fill<-brewer.pal(9,"Set1")[1:3]
newpos <- pos[!duplicated(pos$position),]
newpos$aachange <- newpos$mutation
rownames(newpos) <-newpos$position
duppos <- names(table(pos$position)[table(pos$position)>1])
newpos[duppos,]$aachange <- newpos[duppos,]$aa

sample.gr <- GRanges("chr1", IRanges(newpos$position, width=1, names=newpos$aachange))
sample.gr$label.parameter.rot <- 45
sample.gr$label <- as.character(table(pos$position))
sample.gr$label.col <- "white"
sample.gr$color <-sample.gr$border <- brewer.pal(9,"Set1")[4]
sample.gr$score<- log2(table(pos$position))
labs<-sample.gr$label
labs[labs=="1"]<-""
sample.gr$label <- labs
torb <- sample(c("top", "bottom"), length(sample.gr), replace=TRUE)
sample.gr$SNPsideID <- torb
sample.gr$color[torb=="bottom"] <- sample.gr$border[torb=="bottom"] <- brewer.pal(9,"Set1")[7]

pdf("lolliplot.pdf",width=21)
lolliplot(sample.gr, features,xaxis = c(protein_inform$Start,protein_inform$End),yaxis = F, ylab = F,type="circle")
dev.off()
#Plotting code of Figure4
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
pdata_melt <- read.csv("TMB.csv", header = T)
c <- ggplot(pdata_melt,
            aes(x=Sianature, y=TMB, fill=Type)) +
  geom_boxplot(notch = F, alpha=0.95, outlier.colour = "black",
               outlier.shape = 16,
               outlier.size = 2) +
  scale_fill_manual(values= c("#E7B800","#2E9FDF")) +
  scale_color_manual(values= c("#E7B800","#2E9FDF")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 13), 
        axis.text.y = element_text(angle = 90, size = 13),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")
p <- c + stat_compare_means(label = "p.signif")
ggsave(p, filename = "TMB.pdf", width = 5, height = 7)
#Plotting code of NAL、MANTIS score is consistent with TMB
#Plotting code of Figure5A and Figure5B
#CIBERSORT algorithm
CoreAlg <- function(X, y){

  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)

  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#' do permutations
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #print(itor)

    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

#' Main functions
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mixture_file heterogenous mixed expression
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @export
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  library(e1071)
  library(parallel)
  library(preprocessCore)

  #read in data
  X <- read.table(sig_matrix,header=T,sep="\t",row.names=1,check.names=F)
  Y <- read.table(mixture_file, header=T, sep="\t", row.names=1,check.names=F)

  X <- data.matrix(X)
  Y <- data.matrix(Y)

  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  while(itor <= mixtures){

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1

  }

  #save results
  write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=F, col.names=F, quote=F)

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}
library("limma")  
expFile="symbol.txt"
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
v <-voom(data, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.txt",sep="\t",quote=F,col.names=F)        
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.txt", perm=1000, QN=TRUE)
#Plotting code of Figure5C
d1 <- STADR
for(i in 1:nrow(d1)){
    d1[i,1:i] <- STADR[i,1:i]
}
d2 <- STADP
for (i in 1:nrow(STADP)) {
    STADP[i,1:i] <- STADP[i,1:i]
}
d1[d2 > 0.05] <- NA
library(ComplexHeatmap)
colCorRight <-  circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "#ef3b2c"))
colCorLeft <- circlize::colorRamp2(c(-1, 0, 1), c("yellow", "white", "#762a83"))
p1 <- Heatmap(d1, rect_gp = gpar(type = "none"), 
              show_heatmap_legend = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height,
                          gp = gpar(col = "grey", fill = NA))
                if(i == j) {
                  grid.circle(x = x, y = y, r = 0.5 * min(unit.c(width, height)), gp = gpar(fill = "grey", col = NA))
                  }else if(i > j) {
                    grid.circle(x = x, y = y, r = abs(datR[i, j])/2 * min(unit.c(width, height)), 
                                gp = gpar(fill = colCorLeft(datR[i, j]), col = NA))
                    } else {
                      grid.circle(x = x, y = y, r = abs(datR[i, j])/2 * min(unit.c(width, height)), 
                                  gp = gpar(fill = colCorRight(datR[i, j]), col = NA))
                      }
                },
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = T, show_column_names = T, 
              row_names_side = "right", 
              row_names_rot = 45,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8)
              )
lgdRight <- Legend(col_fun = colCorRight,
                   direction = "horizontal")
lgdLeft <- Legend(col_fun = colCorLeft,
                  direction = "horizontal")
pd = list(lgdRight, lgdLeft) 
pdf("Figure5C.pdf", width = 5, height = 5.5)
draw(p1, annotation_legend_list = pd,
     annotation_legend_side = "top")
dev.off()
#Plotting code of Figure5D and Figure5E
mygene_data <- read.csv("input_expr.csv", row.names = 1)
Subtype <- read.csv("input_group.csv", row.names = 1)
colnames(mygene_data)<-mygene_data_test
com_sam <- intersect(colnames(mygene_data),rownames(Subtype))
mygene_data <- mygene_data[,com_sam]
comprTab <- cross_subtype_compr(expr = mygene_data, 
                                subt = Subtype,
                                two_sam_compr_method = "wilcox", 
                                res.path = ".")
n.show_top_gene <- nrow(mygene_data)
subt.order <- Subtype[order(Subtype$Subtype),,drop = F]
indata <- mygene_data[comprTab$gene[1:n.show_top_gene],rownames(subt.order)]
plotdata <- t(scale(t(indata)))
plotdata[plotdata > 2] <- 2
plotdata[plotdata < -2] <- -2
blank <- "    "
p.value <- comprTab$adjusted.p.value[1:n.show_top_gene]
sig.label <- ifelse(p.value < 0.001,"****",
                    ifelse(p.value < 0.005,"***",
                           ifelse(p.value < 0.01,"**",
                                  ifelse(p.value < 0.05,"*",""))))
p.label <- formatC(p.value, 
                   format = "e",
                   digits = 2) 
add.label <- str_pad(paste0(rownames(plotdata),sig.label), 
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")

annCol <- subt.order 
colnames(annCol)[1] <- paste(str_pad(colnames(annCol)[1], 
                                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                                     side = "right"),
                             "P-value",
                             sep = blank)

annColors <- list(c("WT"="lightblue", "MT"="pink")) 
names(annColors) <- colnames(annCol)[1] 
pheatmap(cellheight = 15, cellwidth = 1,
         mat = plotdata, 
         scale = "none", 
         annotation_col = annCol,
         annotation_colors = annColors, 
         cluster_cols = F,
         cluster_rows = F, 
         show_colnames = F,
         show_rownames = T,
         #annotation_legend = F, 
         labels_row = paste(add.label, p.label, sep=blank), 
         fontfamily = "mono", 
         gaps_col = c(21),
         filename = "all_ICI_heatmapPvalue.pdf")
#Plotting code of  Figure5F-Figure5L
p <- wilcox.test(d1[which(d1$MUC16 == "A"),"Th2.Cells"],d1[which(d1$MUC16 == "B"),"Th2.Cells"])$p.value
jco <- c("#4876FF","#FFA54F")
ggplot(data = d1,aes(x = MUC16, y = Th2.Cells, fill = MUC16))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),
              size=0.8, color="black") + 
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7)+ 
  geom_point(shape = 21, size=2, 
             position = position_jitterdodge(), 
             color="black", alpha=1)+ 
  theme_classic() +
  ylab(expression) +
  xlab("")  +
  annotate(geom="text", cex=6,
           x=1.5, y=1, 
           label=paste0("P ", ifelse(p<0.001, "< 0.001", paste0("= ",round(p,3)))), 
           color="black") + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), 
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))
ggsave("Th2_violin.pdf", width = 4.5, height = 4)
#The plotting codes for other immune scores are the same as above
#Plotting code of Figure6A and Figure6B
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
data <- read.csv("sequence.csv",header=TRUE)
gene <- data$SYMBOL
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
dim(data_all)
head(data_all)
data_all_sort <- data_all %>% 
  arrange(desc(logFC))
head(data_all_sort)
geneList = data_all_sort$logFC 
names(geneList) <- data_all_sort$ENTREZID 
kegg_gmt <- read.gmt("c2.cp.kegg.v7.4.entrez.gmt")
reactome_gmt<-read.gmt("c2.cp.reactome.v7.4.entrez.gmt")
gsea <- GSEA(geneList,
             TERM2GENE = kegg_gmt)
 gse.GO <- gseGO(
  geneList, 
  ont = "ALL",  
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)
#Visualization via CAMOIP
#Plotting code of Figure6C
library(pRRophetic)
