

if(T){
  library(data.table)
  library(dplyr)
  library(survival)
  library(survminer)
  library(stringr)
  library(ggplot2)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
}


#### unicox and multi cox ####
load("D:/Aproj/aBULK/TCGA_LUAD_dealed0104.rda")

sur <- cl2[,c('id','os','os_time')]
sur <- na.omit(sur)

mat <- as.data.frame(t(zscore))
mat <- mat[,sur$id]


dat = mat[setdiff(rownames(mat),grep("-", rownames(mat), value = TRUE)),] 

surexp <- cbind(sur,as.data.frame(t(dat)))

unicox = function(y){
  covariates<-colnames(y)[3:ncol(y)] 
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(os_time,os)~', x)))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = y)}) 
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                          
                           p.value<-signif(x$wald["pvalue"], digits=3)
                           
                           HR <-signif(x$coef[2], digits=3);
                         
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           Hazard_Ratio <- paste0(HR, " (", 
                                                  HR.confint.lower,"-",HR.confint.upper, ")")
                           res<-c(p.value,HR,HR.confint.lower,HR.confint.upper,Hazard_Ratio)  
                           names(res)<-c("p.value","HR","lower","upper","HR (95% CI for HR)")
                           return(res)
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  result=as.data.frame(res)
  result= na.omit(result)
  result = result %>% arrange(p.value) #%>% 
    #filter(p.value<0.05)
  result$HR=as.numeric(result$HR)
  result =result %>% arrange(HR)
  res= result
  return(res)
}

unires = unicox(surexp)


write.csv(cl2,file="./out/TCGA_cl.csv")

clnew <- data.table::fread("D:/Aproj/aBULK/out/TCGA_cl.csv",data.table = F) 

table(clnew$stage=='NA')

multi <- clnew[match(sur$id,clnew$id),]
identical(multi$id,sur$id)

muldat <- cbind(multi[,c('id','os','os_time','stage','age','gender')],as.data.frame(t(dat)) )

muldat$gender <- factor(muldat$gender)
muldat$stage <-  dplyr::case_when(
  muldat$stage %in% c('Stage I','Stage IA','Stage IB') ~ 'Stage I',
  muldat$stage %in% c('Stage II','Stage IIA','Stage IIB')~ "Stage II",
  muldat$stage %in% c('Stage IIIA','Stage IIIB') ~ "Stage III",
  muldat$stage %in% c('Stage IV') ~ "Stage IV"
                                   )

muldat$stage <- factor(muldat$stage)

multicox = function(y){
  covariates<-colnames(y)[7:ncol(y)] 
  covariates <-'DNER'
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(os_time,os)~', x,'+age + gender + stage')))
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = muldat)}) 
  
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x1 <- summary(x)
                       
                           p.value<-signif(x1$coefficients[1,'Pr(>|z|)'], digits=3)
                          
                           HR <-signif(x1$conf.int[1,'exp(coef)'], digits=3);
                          
                           HR.confint.lower <- signif(as.numeric(x1$conf.int[1,'lower .95']), 3)
                           HR.confint.upper <- signif(as.numeric(x1$conf.int[1,'upper .95']), 3)
                           
                           Hazard_Ratio <- paste0(HR, " (",as.character(HR.confint.lower),"-",as.character(HR.confint.upper), ")")
                           res<-c(p.value,HR,HR.confint.lower,HR.confint.upper,Hazard_Ratio) 
                           names(res)<-c("p.value","HR","lower","upper","HR (95% CI for HR)")
                           return(res)
                         })
  
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  result=as.data.frame(res)
  result= na.omit(result)
  result = result %>% arrange(p.value) #%>% 
  #filter(p.value<0.05)
  result$HR=as.numeric(result$HR)
  result =result %>% arrange(HR)
  res= result
  return(res)
}



mulres <- multicox(muldat)

univ_models
ggforest(univ_models$DNER,
         data= muldat,
         main= "Hazard ratio",
         cpositions= c(0.02, 0.22, 0.4), 
         fontsize= 1, 
         refLabel= "1", 
         noDigits= 3)+  
scale_fill_manual(values = c("#FFFFFF33",
                             "#00000033"), guide = "none")

library(forestmodel)
forest_model(univ_models$DNER)
devtools::install_github("kassambara/coxplot")
library(coxplot)
coxforest(fit)


unires005 <- unires %>% filter(p.value<0.05) %>% mutate(gene=rownames(unires005))
mulres005 <- mulres %>% filter(p.value<0.05) %>% mutate(gene=rownames(mulres005))

table(duplicated(unires005$gene))
table(duplicated(mulres005$gene))

rownames(unires005) <- NULL
rownames(mulres005) <- NULL

inter <- intersect(unires005$gene,mulres005$gene)

n1 <- unires005[(unires005$gene%in%inter),] %>% arrange(gene)

n2 <- mulres005[(mulres005$gene%in%inter),] %>% arrange(gene)

uni_mul <- cbind(n1,n2) 

colnames(uni_mul)[7:12] <- paste0(colnames(uni_mul)[7:12],'_mul')

uni_mul <- uni_mul %>% arrange(desc(HR))

#### 
save(unires,unires005,mulres,
     mulres005,uni_mul,muldat,surexp,file="./out/TCGA_uni_multi_reg_res.rda")


uni_mul_resTCGA <- list(
  unires = unires,
  unires005 = unires005,
  mulres = mulres,
  mulres005 = mulres005,
  uni_mul = uni_mul,
  unidat = surexp,
  muldat = muldat
)

save(uni_mul_resTCGA, file = "./out/TCGA_uni_multi_reg_res.rda")


#### ggplot2  ####
library(ggplot2)

mat_long <- reshape2::melt(top_anno2[,c(1,7:53)])

mytheme <-theme(
  panel.grid = element_blank(), #去网格
  panel.background = element_rect(fill = "transparent",colour = 'white',linewidth=1),#空白背景 fill 背景填充
  plot.title = element_text(lineheight=1,size=24),
  plot.subtitle = element_text(lineheight=1,size=21),
  # axis.line = element_line(arrow = arrow(length = unit(0.3, "cm")),linewidth=0.8),#箭头
  # axis.line.x = element_line(colour = "black") ,
  # axis.line.y = element_line(colour = "black") ,
  # axis.text.x = element_text(angle=45,hjust = 1,size=16),
  axis.text.y = element_blank(),
  axis.text.x = element_text(size=14,angle=60,hjust=0,family = "sans"),
  axis.title.x = element_text(angle=0,size=14),
  axis.title.y = element_text(angle=90,size=14),
  #axis.ticks = element_blank, 
  axis.ticks = element_blank(),
  #axis.ticks = element_line(colour = "black",linewidth = 1),
  axis.ticks.length.x = unit(4, "pt"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14,family = "sans"), #sans Arial
  #legend.position = "none",
  legend.position =  "right",
  legend.key.size = unit(2, "lines"), 
  
  legend.key = element_rect(fill = "white",colour = 'black',linewidth=1.2) #### 图例元素边框
)

# levels(mat_long$variable)
# table(mat_long$variable)

# 
ggplot(mat_long, aes(y = id, x = variable, fill = value)) +
  geom_tile(linejoin='round',alpha=1,width=2) +
  mytheme + scale_x_discrete(position = "top")+
  scale_fill_viridis_b(direction = 1)+labs(x='Risk signatures')



if(T){
  load("D:/Aproj/aBULK/out/TCGA_uni_multi_reg_res.rda")
  load("D:/Aproj/aBULK/out/GSE13213_uni_multi_reg_res.rda")
  load("D:/Aproj/aBULK/out/GSE31210_uni_multi_reg_res.rda")
  load("D:/Aproj/aBULK/out/GSE68465_uni_multi_reg_res.rda")
  load("D:/Aproj/aBULK/out/GSE50081_uni_multi_reg_res.rda")
}

dat_list <- list( uni_mul_res13213=uni_mul_res13213,
                  uni_mul_res31210=uni_mul_res31210,
                  uni_mul_res50081=uni_mul_res50081, 
                  uni_mul_res68465=uni_mul_res68465,
                  uni_mul_resTCGA=uni_mul_resTCGA
                  )





#### KM plot ####

if(T){
  library(survival)
  library(survminer)
  library(caret)
  library(survex)
  library(ranger)
}

mythe = theme_bw()+theme(
  #panel.grid = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = 'black'),
  axis.line.x = element_line(colour = "black") , #轴线宽
  axis.line.y = element_line(colour = "black") ,
  axis.text.x = element_text(colour = "black",size = 16), #轴标注
  axis.text.y = element_text(colour = "black",size = 16),
  axis.ticks = element_line(color = "black", size = 1), 
  axis.ticks.length = unit(5, "pt"), 
  axis.title.x = element_text(colour = "black",size = 16),
  axis.title.y = element_text(colour = "black",size = 16),
  plot.title = element_text( size = 16),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14,family = "sans")
  
  #legend.position =  "right",
  
)

y = 'sss'

psur= function(x,y){
  cut <- surv_cutpoint(uni_mul_res72094$unidat,'os_time','os',y)
  cut
  ## 
  cat <- surv_categorize(cut)
  fitsur <- survfit(Surv(os_time,os)~EFNA3,cat)
  n1=nrow(subset(cat,EFNA3=="high"))
  n2=nrow(subset(cat,EFNA3=="low"))
  
  ggsurvplot(
    fitsur, cat,
    palette = c('purple','yellow3'),
    size = 1.5,
    censor = TRUE,
    censor.shape = "|", censor.size = 6,
    pval = F ,
    linetype = 1,
    conf.int = TRUE, 
    conf.int.style = 'step',
    break.x.by = 1000,
    break.y.by = 0.25,
    pval.method = FALSE,
    pval.size = 10,
    risk.table = FALSE,
    legend = "bottom",
    legend.labs = c(paste0("High \n(n=", n1, ")"), paste0(" Low \n(n=", n2, ")")),
    legend.title = 'EFNA3',
    xlab = "Days",
    ylab = 'OS',
    ggtheme = theme_bw() + 
      theme(legend.position = "bottom",
            legend.background = element_rect(fill="white", colour="black", size=0.5, linetype="solid"),
            plot.title = element_text(size=18),
            axis.title = element_text(size=14),
            axis.text = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=13))
  ) + ggtitle("Log-rank: p < 0.0001")
          
  
}



options(warning=-1)
if(T){
  library(tidyverse)
  library(arrow)
  library(data.table)
  library(fastmatch)
  library(coloc)
  library(TwoSampleMR)
  library(gwasvcf)
  library(gwasglue)
  library(BSgenome)
  library(MungeSumstats)
  library(VariantAnnotation)
  library(LDlinkR)
  #install.packages('D:/R432/Rbase432/R-4.3.2/library/SNPlocs.Hsapiens.dbSNP155.GRCh38_0.99.24.tar.gz')
}


#### CNV  RACGAP1    ####
if(T){
  library(survival)
  library(survminer)
  library(caret)
  library(survex)
  library(ranger)
  library(tidyverse)
}

if(T){
  load("D:/Aproj/aBULK/out/TCGA_uni_multi_reg_res.rda")
  load("D:/Aproj/aBULK/Data/GSVA/LUAD_cnv_sp516.rda")
}

mrna <- uni_mul_resTCGA$unidat
cnv <- as.data.frame(t(cnv))
cnv2 <- cnv %>% mutate(id = rownames(cnv)) %>% select(id,everything()) %>% 
  select(c('id','RACGAP1')) %>% 
  rename(racgap1='RACGAP1')

mer <- mrna %>% left_join(cnv2,by='id')

p <- ggplot(mer, aes(x = racgap1, y = RACGAP1)) +
  geom_jitter(shape=16,alpha = 1,width = .3,
              #fill='pink',
              color = 'green4', size = 2) +
  geom_smooth(method = "lm", size=2,se = TRUE, color = 'blue', formula = y ~ x) +
  labs(
    title = "Correlation between RACGAP1 CNV and RACGAP1 mRNA",
    x = "CNV RACGAP1",
    y = "RACGAP1 mRNA",
    caption = "rho=0.40, p<0.0001"
  ) +
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
p
ggsave(width=3,height=3,file='D:/Aproj/aMR/Pt_lung1107/cnv_luad/Bulk_tcga_mrna_RACGAP1_cnv_cor_exp.pdf')


library(ggpubr)

cat <- mer[,c('os_time','os','racgap1')] %>% 
  rename(cnv='racgap1') %>% 
  filter(cnv %in%c(1,-1) ) #%>% 
#filter(os_time<=2000)

fitsur <- survfit(Surv(os_time, os) ~ cnv, data = cat)

p <- ggsurvplot(
  fitsur, data = cat,
  pval = T,
  size = 1,
  palette = c('blue', 'tomato'),
  censor = TRUE,
  censor.shape = "|", censor.size = 1,
  conf.int = FALSE, 
  risk.table = T,
  legend = "bottom",
  xlab = "Days",
  ylab = 'Percent Survival',
  ggtheme = theme_pubr(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      legend.background = element_rect(fill="white", color="black", size=0.5),
      plot.title = element_text(size=12),
      axis.title = element_text(size=12),
      axis.text = element_text(size=10),
      legend.text = element_text(size=12),
      legend.title = element_text(size=12)
    )
)

p_combined <- ggarrange(p$plot, p$table, nrow = 2, heights = c(2, 1))
p_combined 
ggsave(
  'D:/Aproj/aMR/Pt_lung1107/cnv_luad/Bulk_tcga_mrna_RACGAP1_cnv_sur_delvsplus_final.pdf',
  plot = p_combined, 
  width = 3.2,
  height = 4
)



#### survival ####

if(T){
  library(survival)
  library(survminer)
  library(caret)
  library(survex)
  library(ranger)
}


var_name <-gene1[i]  
cat <- uni_cox_exp[,c('time','event')]        
cat[[var_name]] <- ifelse(uni_cox_exp[[var_name]] > median(uni_cox_exp[[var_name]]), 'high', 'low')
cat[[var_name]] <- factor(cat[[var_name]], levels = c("high", "low"))

n1 <- sum(cat[[var_name]] == "high")
n2 <- sum(cat[[var_name]] == "low")

fitsur <- survfit(Surv(time, event) ~ get(var_name), data = cat)

p <- ggsurvplot(
  fitsur, data = cat,
  legend.title = var_name,
  palette = c('#785EF0', '#D0D0D0'),
  pval = F,
  size = 1.2,
  censor = TRUE,
  censor.shape = "|", censor.size = 2.5,
  conf.int = FALSE, 
  risk.table = FALSE,
  legend = "bottom",
  legend.labs = c(paste0("High (n=", n1, ")"), paste0("Low (n=", n2, ")")),
  xlab = "Days",
  ylab = 'Percent Survival',
  ggtheme = theme_pubr(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      legend.background = element_rect(fill="white", color="black", size=0.5),
      plot.title = element_text(size=12),
      axis.title = element_text(size=12),
      axis.text = element_text(size=10),
      legend.text = element_text(size=12),
      legend.title = element_text(size=12)
    )
)
p$plot <- p$plot + labs(
  #title = paste0(var_name,'\nhigh(n =279),low(n=280)'),
  title  = var_name,
  subtitle = paste0('Logrank p = ', signif(surv_pvalue(fitsur, data=cat)$pval, 2))
)
#print(p)

plist[[i]] <- p$plot



ggpubr::ggarrange(plist[[1]],plist[[2]],plist[[3]],plist[[4]],
                  plist[[5]],plist[[6]],plist[[7]],plist[[8]],
                  plist[[9]],plist[[10]],plist[[11]],plist[[12]],
                  plist[[13]],plist[[14]],plist[[15]],plist[[16]],
                  plist[[17]],plist[[18]],ncol=6,nrow=3,
                  common.legend = T
                  
)
gene1
ggpubr::ggarrange(plist[[1]],plist[[2]],plist[[3]],plist[[4]],
                  plist[[5]],plist[[6]],plist[[7]],plist[[8]],
                  plist[[9]],plist[[10]],plist[[11]],plist[[12]],
                  plist[[13]],plist[[14]],plist[[15]],plist[[16]],
                  plist[[17]],plist[[18]],ncol=6,nrow=3,
                  common.legend = T
                  
)




####  protein RNA compare  ####
sur <- data.table::fread(data.table = F,'D:/Aproj/aMR/Pt_lung1107/TCGA_data_gsca/ucscxena/TCGA-LUAD.survival.tsv')
expr <- data.table::fread(data.table = F,"E:/LUAD_Data/TCGA-LUAD.htseq_fpkm-uq.tsv.gz")
rownames(expr) <- expr$Ensembl_ID
expr <- expr[,-1]

mrna <- trans_exp(expr,mrna_only = T)

gp <- make_tcga_group(mrna)

mrna2 <- sam_filter(mrna)
cl_mat <- match_exp_cl(mrna2,sur,'_PATIENT')
meta <- cl_mat$cl_matched
colnames(meta) <- c("id","sample_id","sample","event","time" )
exp <- cl_mat$exp_matched[,]
intersect(mr_fdr$gene,rownames(exp))
exp2 <- cl_mat$exp_matched[mr_fdr$gene,]

uni_cox_exp <- cbind(meta,t(exp2))

uni_cox_exp[['status_gp']] <- ifelse( substr(uni_cox_exp$sample_id,14,15) %in%'01','tumor','ctl')

library(tidyplots)
library(gghalves)

pdat <- uni_cox_exp

var <- mr_fdr$gene[1]


pbox <- function(var){
  p <- pdat %>% 
    filter(!!sym(var) !=0 ) %>% 
    ggplot( aes(x = status_gp, y = !!sym(var), fill = status_gp)) +
    geom_half_boxplot(alpha=.5,color='black',size=.8,width = 1) +
    #geom_jitter(aes(x=gp,y=!!sym(var)),alpha = 1,size=2)+
    theme_bw(base_size = 16)+
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          axis.text.x = element_text(size=16,angle=45,hjust=1),
          plot.subtitle = element_text(face='italic')
    )+
    scale_fill_manual(values=c('purple','yellow') )+
    labs(subtitle = var,x='',y='fpkm') +
    ggpubr::stat_compare_means(
      comparisons = list(c("tumor", "ctl")),
      method = "wilcox.test",
      label = "p.signif",                 # 显示格式化的 p 值
      bracket.size = .5,
      size=7.5,
      vjust = 1.5,
      tip.length = 0.0
    )
  p
  
}

gene <- mr_fdr$gene

plist <- lapply(gene, pbox)

plist2 <- ggpubr::ggarrange(
  plist[[1]],plist[[2]],plist[[3]],plist[[4]],plist[[5]],
  plist[[6]],plist[[7]],plist[[8]],plist[[9]],plist[[10]],
  plist[[11]],plist[[12]],plist[[13]],plist[[14]],
  plist[[15]],plist[[16]],plist[[17]],plist[[18]],nrow=2,ncol=9
)

#plist2 
ggsave(plist2,width = 14,height = 6,
       file= paste0( "D:/Aproj/aMR/Pt_lung1107/output_new/tcga_exp_tumor_vs_ctl/",'mer','new.pdf'))



#### vad ####
if(T){
  
  load("D:/Aproj/aBULK/out/GSE13213_uni_multi_reg_res.rda")
  # load("D:/Aproj/aBULK/out/GSE31210_uni_multi_reg_res.rda")
  # load("D:/Aproj/aBULK/out/GSE68465_uni_multi_reg_res.rda")
  # load("D:/Aproj/aBULK/out/GSE50081_uni_multi_reg_res.rda")
}

inter <- intersect(gene,rownames(uni_mul_res13213$unires))
uni_mul_res13213$unires005['RACGAP1',]


slc <- c("RACGAP1","CDH17","CTSH","GOLM1")



get__df <- function(model_result, genes) {
  
  slc <- c("RACGAP1","CDH17","CTSH","GOLM1")
  uni <- model_result$unires[slc, ]
  uni[["id"]] <- paste0(rownames(uni), "_uni")
  mul <- model_result$mulres[slc, ]
  mul[["id"]] <- paste0(rownames(mul), "_mul")
  
  uni_mul <- rbind(uni, mul)
  uni_mul <- uni_mul %>% arrange(id)
  
  return(uni_mul)
}

g32 <- get__df(uni_mul_res31210)

df <- g32 %>%
  mutate(across(c(HR, lower, upper, p.value), as.numeric))


df$Label <- df$id
df$Label <- factor(df$Label, levels = rev(df$Label))
df$Label2 <- paste0(df$Label, "    p=", sprintf("%.3g", df$p.value))

p <- ggplot(df, aes(x = HR, y = Label2)) +
  
  geom_segment(aes(x=lower, xend=upper, y=Label2, yend=Label2), color="#377eb8", size=2) +
  geom_point(color="black", fill="white", shape=21, size=5, stroke=1.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +

  scale_x_continuous(trans = "log10", limits = c(0.4, 2.5), 
                     breaks = c(0.5, 0.7, 1, 1.5, 2)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size=12, color="black"),
    axis.title.x = element_text(size=13, margin = margin(t=10)),
    plot.title = element_text(size=15, face="bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85"),
    axis.line.x.bottom = element_line(color="black"),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x = "Hazard Ratio (log scale)",
    title = "Forest plot of Hazard Ratios"
  )
p

ggsave(p,width = 5,height = 3,
       file= paste0( "D:/Aproj/aMR/Pt_lung1107/output_new/",'gse31210','f.pdf'))


#### 
g50 <- get__df(uni_mul_res50081)

df <- g50 %>%
  mutate(across(c(HR, lower, upper, p.value), as.numeric))


df$Label <- df$id
df$Label <- factor(df$Label, levels = rev(df$Label))
df$Label2 <- paste0(df$Label, "    p=", sprintf("%.3g", df$p.value))

p <- ggplot(df, aes(x = HR, y = Label2)) +
  
  geom_segment(aes(x=lower, xend=upper, y=Label2, yend=Label2), color="#377eb8", size=2) +
  geom_point(color="black", fill="white", shape=21, size=5, stroke=1.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  # p值显著性高亮
  # geom_text(aes(x = 1, label = paste0("p=", sprintf("%.3g", p.value))), 
  #           hjust = 0, 
  #           color = ifelse(df$p.value<0.05, "black", "black"),
  #           size = 4) +
  scale_x_continuous(trans = "log10", limits = c(0.4, 2.5), 
                     breaks = c(0.5, 0.7, 1, 1.5, 2)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size=12, color="black"),
    axis.title.x = element_text(size=13, margin = margin(t=10)),
    plot.title = element_text(size=15, face="bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85"),
    axis.line.x.bottom = element_line(color="black"),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x = "Hazard Ratio (log scale)",
    title = "Forest plot of Hazard Ratios"
  )
p

ggsave(p,width = 5,height = 3,
       file= paste0( "D:/Aproj/aMR/Pt_lung1107/output_new/",'gse50081','f.pdf'))




if(T){
  library(survival)
  library(survminer)
  library(corrplot)
  library(glmnet)
  library(caret)
  library(leaps)
  library(tidyverse)
  library(randomForestSRC)
  library(ggplot2)
  library(survex)
  library(ranger)
  library(ggpubr)
}


#### CPTAC diff  ####

library(tidyverse)
library(tidyplots)
library(rlang)
library(ggsci)

slc <- c("RACGAP1","CDH17","CTSH","GOLM1")

if(T){
  tum <-arrow::read_tsv_arrow("D:/Aproj/aMR/Pt_lung1107/cptac/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt")
  ctl <- arrow::read_tsv_arrow("D:/Aproj/aMR/Pt_lung1107/cptac/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt")
  tum$idx <- substr(tum$idx,1,15)
  ctl$idx <- substr(ctl$idx,1,15)
  tum <- tum[!duplicated(tum$idx),]
  ctl <- ctl[!duplicated(ctl$idx),]
  inter <- intersect(ctl$idx,tum$idx)
  ctl <- ctl[ctl$idx %in%inter ,]
  tum <- tum[tum$idx %in%inter ,]
  names(tum)[-1] <- paste0(names(tum)[-1],'_tum')
  dat <- cbind(tum,ctl)
  dat[is.na(dat)] <- 0
  
  load("D:/Anno files and signatures/Gene_data_for_annotion.Rdata")
  ID$gene_id <- substr(ID$gene_id,1,15)
  ID <- ID[!duplicated(ID$gene_id),] |> 
    dplyr::select(c('gene_id','gene_name'))
  names(dat)[1] <- 'gene_id'
  pt <- ID |> dplyr::right_join(dat,by='gene_id')
  pt <- pt |>dplyr::select(-c('idx'))
  pt <- as.data.frame(t(pt))
  pt <- pt[-1,]
  colnames(pt) <- pt[1,]
  pt <- pt[-1,]
  pt$idx <- rownames(pt)
  pt$gp <- ifelse(grepl('tum',rownames(pt)),'Tum','ctl' )
  pt$idx
  pt$gp
  pt[,c(ncol(pt)-1,ncol(pt))]
  
  
}

var <- 'RACGAP1'
var <-  slc[4]
pdat <- pt %>% select(c('gp',!!sym(var)))    
pdat[[var]] <- as.numeric(pdat[[var]])


pdat[['sample']] <- rownames(pdat)

pdat$sample  <- sub('_tum$','',pdat$sample,)

p <- pdat %>% 
  filter(!!sym(var) !=0 ) %>% 
  ggplot( aes(x = gp, y = !!sym(var), fill = gp)) +
  geom_half_violin(alpha=.5,color='black',size=.8,width = 1) +
  #geom_jitter(aes(x=gp,y=!!sym(var)),alpha = 1,size=2)+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c('lightgrey','purple') )+
  labs(x='Group',y='log2(MS1 intensity)') +
  ggpubr::stat_compare_means(
    comparisons = list(c("ctl", "Tum")),
    method = "wilcox.test",
    label = "p.signif",                 # 显示格式化的 p 值
    bracket.size = .5,
    size=7.5,
    vjust = 1.5,
    tip.length = 0.0
  )
p

ggsave(p,width = 2.7,height = 3,
       file= paste0( "D:/Aproj/aMR/Pt_lung1107/prot_cptac",var,'new.pdf'))



var <-  slc[4]

pdat <- pt %>% select(c('gp',!!sym(var)))    
pdat[[var]] <- as.numeric(pdat[[var]])

pdat[['sample']] <- rownames(pdat)
pdat$sample  <- sub('_tum$','',pdat$sample,)

pdat2 <- pdat %>% filter(!!sym(var)!=0)
wide_data <- pdat2 %>%
  select(gp, !!sym(var), sample) %>%
  pivot_wider(names_from = gp, values_from = !!sym(var)) %>% 
  na.omit()

t_result <- t.test(wide_data$ctl, wide_data$Tum, paired = TRUE)

pvalue <- sprintf("%.3e",t_result$p.value )

p <- ggplot(pdat2) +
  
  geom_line(aes(x = gp, y = !!sym(var),fill=gp, group = sample),color = "gray60", size = 0.7, alpha = 0.7)+
  geom_jitter( aes(x = gp, y = !!sym(var),color=gp, group = sample),
               shape=21,size=1,width =0,alpha=.7) +
  scale_color_manual(values=c("purple", "pink4"))+
  geom_boxplot( size=1,aes(x = gp, y = !!sym(var)),fill=NA,width=.35 )+
  #scale_color_nejm()+
  labs(x='')+
  theme_classic(base_size = 14)+
  ggtitle(paste0('Paired t-test = ',pvalue))+
  theme(legend.position = 'none')

p

ggsave(p,width = 3,height = 3.2,
       file= paste0( "D:/Aproj/aMR/Pt_lung1107/output_new/CPTAC/",var,'match.pdf'))




sur <- arrow::read_tsv_arrow("D:/Aproj/aMR/Pt_lung1107/cptac/LUAD_survival.txt")
meta <- arrow::read_tsv_arrow("D:/Aproj/aMR/Pt_lung1107/cptac/LUAD_meta.txt")





#### bulk mRNA  RA corr  ####

load("D:/Aproj/aBULK/TCGA_LUAD_dealed0104.rda")

results <- apply(zscore, 2, function(x) {
  test_result <- cor.test(as.numeric(zscore[,'RACGAP1']),as.numeric(x), method = 'spearman')
  return(c(correlation = test_result$estimate, p_value = test_result$p.value))
})
results_df <- as.data.frame(t(results), stringsAsFactors = FALSE)
names(results_df) <- c("Correlation", "P_Value")
results_df <- results_df %>% arrange(Correlation) %>% 
  mutate(FDR=p.adjust(P_Value))

df <- results_df %>%  arrange(desc(Correlation))


load(file="D:/Aproj/aMR/Pt_lung1107/cnv_luad/Tumor_vs_normal__ra_plus_vs_del_tumor.rda")

diff_genes <- markers2 %>% filter(pt_diff>0.1&p_val_adj<0.05)
cor_genes <- df %>% filter(Correlation>0.5,FDR<0.05)

inter <- intersect( rownames(cor_genes),rownames(diff_genes) ) 

write.csv(inter,
          file='D:/Aproj/aMR/Pt_lung1107/PPI_core_prot_gene.csv')

####  GSEA  correlation cor ####

if(T){
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(tidyverse)
  library(GseaVis)
  library(msigdbr)
}

deff <- df %>% dplyr::select(c('Correlation','P_Value')) %>% 
  rename(logFC='Correlation',P.Value='P_Value')

if(T){
  reslimma = data.frame(gene =rownames(deff ),logFC=deff$logFC)
  ranklist = function(x){
    x= reslimma
    entrez <- bitr(x$gene,
                   fromType = "SYMBOL",#
                   toType = "ENTREZID",#型
                   OrgDb = "org.Hs.eg.db")
    
    entrez <- entrez[!duplicated(entrez$SYMBOL),]
    list = x[x$gene%in%entrez$SYMBOL,]
    if(identical(list$gene,entrez$SYMBOL)){
      print("construct gene list sucessfully")
    }
    list1 = cbind(list,entrez)
    if(identical(list1$gene,list1$SYMBOL)){
      print("construct merge list sucessfully")
    }
    
    list1 = list1[order(list$logFC,decreasing = T),]
    generank = as.vector(as.numeric(list1$logFC))
    names(generank) = list1$ENTREZID
    return(generank)
  }
  genelist = ranklist(reslimma)
  ranklist2 = function(x){
    x= reslimma
    entrez <- bitr(x$gene,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = "org.Hs.eg.db")
    
    entrez <- entrez[!duplicated(entrez$SYMBOL),]
    list = x[x$gene%in%entrez$SYMBOL,]
    if(identical(list$gene,entrez$SYMBOL)){
      print("construct gene list sucessfully")
    }
    list1 = cbind(list,entrez)
    if(identical(list1$gene,list1$SYMBOL)){
      print("construct merge list sucessfully")
    }
    
    list1 = list1[order(list$logFC,decreasing = T),]
    generank = as.vector(as.numeric(list1$logFC))
    names(generank) = list1$SYMBOL
    return(generank)
  }
  list3 =  ranklist2(reslimma)
}

m_df = msigdbr(species = 'Homo sapiens' , category = "H") %>%
  dplyr::select(gs_name,gene_symbol) 

colnames(m_df)<-c("term","gene")


h_ges_bulkra <- GSEA(list3, 
                     TERM2GENE = m_df,
                     minGSSize = 5,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     verbose = FALSE,
                     eps = 0)


id <- h_ges_bulkra@result %>% rownames()

id

p <- gseaNb(object = h_ges_bulkra,
            subPlot = 2,
            base_size = 14,
            addPval=T,
            geneSetID = id[14])

p     

pth <- "D:/Aproj/aMR/Pt_lung1107/output_new/gsea_TCGA_ra_high_cor_vs_ra_low"

ggsave(p,width=3,height=2.2,file=paste0(pth,'_dna_repair.pdf')
)

sigg <- c('MKI67','CDK1','CDC20','CCNB2','CDCA3')
sigg <- c('PRC1','KIF14','KIF23','KIF4A','ANLN','ECT2','GMNN')
sigg <- c('BUB1','BUB1B','BUB3','MAD2L1','TTK','PLK1','AURKA','AURKB')

cin70 <- c("TPX2", "PRC1", "FOXM1", "CDC2", 
           "TGIF2", "MCM2", "H2AFZ",
           "TOP2A", "PCNA", "UBE2C", "MELK", 
           "TRIP13", "CNAP1", "MCM7", "RNASEH2A",
           "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1",
           "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", 
           "RFC4", "ATAD2", "ch-TOG", "NUP205", "CDC20"
           , "CKS2", "RRM2", "ELAVL1", "CCNB1", "RRM1", 
           "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5", 
           "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "BRRN1",
           "MCM10", "LSM4", "MTB", "ASF1B", "ZWINT", "TOPK", 
           "FLJ10036", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2",
           "RAD21", "ACTL6A", "GPIandMGC13096", "SFRS2",
           "HDGF", "NXT1", "NEK2", "DHCR7", "STK6", "NDUFAB1",
           "NEMP1", "KIF4A")


pdt <- df |> mutate(gene=rownames(df)) |>
  filter(gene %in% cin70 )

order <- pdt$gene

pdt$gene <- factor(pdt$gene, levels = order)

p <- ggplot(pdt, aes(x = gene, y = Correlation)) +
  geom_col(width=1, alpha=.5, fill='red4', color='black') +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle=60, hjust=1)) +
  geom_hline(yintercept=0.5,size=1, linetype="dashed", color="black")

p


ggsave(p,width=14,height=3,
       file="D:/Aproj/aMR/Pt_lung1107/output_new/RA相关性/ra_cell_stab_all.pdf")








