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
