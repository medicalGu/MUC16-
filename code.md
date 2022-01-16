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
