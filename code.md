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
#KM plot(Figure2C)
library("survival")
library("survminer")
