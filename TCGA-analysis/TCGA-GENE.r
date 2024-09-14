#TCGA GENE
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(rstatix) 
library(ggprism)
expression <- read.table('D:\\singleDATA\\TCGA\\expression.csv',header = T,fill = T,sep = ",")
load('TCGA-THCA.Rdata')

group_list <- ifelse(as.numeric(str_sub(rownames(gene_exp),14,15)) < 10,'tumor','normal')
group_list <- factor(group_list,levels = c("normal","tumor"))
table(group_list)
group_list2 <- ifelse(as.numeric(str_sub(rownames(gene_exp),14,15)) == 1,'T1',
                     ifelse(as.numeric(str_sub(rownames(gene_exp),14,15)) == 2,'T2',
                            'T3','T4'))
group_list2 <- factor(group_list2,levels = c('T1','T2','T3','T4'))
table(group_list2)
group_list3 <- ifelse(as.numeric(str_sub(rownames(gene_exp),14,15)) < 10,'N0','N1')
group_list3 <- factor(group_list,levels = c("N0","N1"))
table(group_list)

#自定义配色：
mycol <- c("#377EB8","#E41A1C")

#计算样本数：
num_normal <- length(group_list[group_list == 'normal'])
num_cancer <- length(group_list[group_list == 'tumor'])

#绘图：
p1 <- ggplot(data = dt,
             aes(x = group1, y = log(dt[,1] + 1), color = group1, fill = group1)) +
  geom_boxplot(alpha = 0.4) +
  geom_jitter(width = 0.3,alpha = 0.6) +
  theme_classic() +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  labs(
    title = paste0('Expression of ',gene,' in ',cancer,' based on Sample types'),
    subtitle = paste0('normal(n = ',num_normal,'), tumor(n = ',num_cancer,')'),
    x = 'TCGA Groups',
    y = 'Log (TPM + 1)'
  ) +
  theme(legend.title = element_blank())
p1

#自定义配色：
mycol2 <- c("#6388B4","#FFAE34","#EF6F6A")

#计算样本数：
T1<- length(group_list2[group_list2 == 'T1'])
T2<- length(group_list2[group_list2 == 'T2'])
T3<- length(group_list2[group_list2 == 'T3'])
T4<- length(group_list2[group_list2 == 'T4'])

#绘图：
p2 <- ggplot(data = dt,
             aes(x = group2, y = log(dt[,1] + 1), color = group2, fill = group2)) +
  geom_boxplot(alpha = 0.4) +
  geom_jitter(width = 0.3, alpha = 0.6) +
  theme_classic() +
  scale_color_manual(values = mycol2) +
  scale_fill_manual(values = mycol2) +
  labs(
    title = paste0('Expression of ',gene,' in ',cancer,' based on Sample types'),
    subtitle = paste0('T1(n = ',T1,'), T2(n = ',T2,'), T3(n = ',T3,'), T4(n = ',T4,')'),
    x = 'TCGA Groups',
    y = 'Log (TPM + 1)'
  ) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = -20, hjust = 0)) 
p2

#t检验(wilcox方法一致)：
dt_t1 <- dt %>% t_test(formula = TPM ~ group1,
                       comparisons = list(c("tumor", "normal")), #指定比较组
                       paired = F) %>% add_significance()

dt_t2 <- dt %>% t_test(formula = TPM ~ group2,
                       comparisons = list(c('T1','T2'),c('T1','T3'),c('T1','T4'),c('T2','T3'),c('T2','T4'),c('T3','T4')),
                       p.adjust.method = "bonferroni",
                       paired = F) %>% add_significance()
dt_t1;dt_t2

p1 + add_pvalue(dt_t1,
                bracket.size = 0.5,
                label.size = 4,
                color = 'black',
                y.position = max(log(dt$TPM)) + 0.3)