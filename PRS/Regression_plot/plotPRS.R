
library(tidyverse)    
library(broom)        
library(writexl)      
library(ggpubr)       
library(viridis)      
library(data.table)
library(dplyr)
plot_list=list()

get_model<- function(prs_path,outcome_threshold){
  sp_score<-prs_path
  
  df_prs<-fread(sp_score) %>% as.data.frame()
  PHENO_FILE = "./pheno.txt"
  df_pheno=fread(PHENO_FILE,header=FALSE) %>% as.data.frame()
  colnames(df_pheno)<-c('IID','FID','outcome')
  df_prs_pheno<- merge(df_prs,df_pheno,by.x=c("IID"),by.y='IID')
  
  prs_outcome <- read.csv("./covar.txt",sep='\t')
  
  df_prs_pheno<- merge(df_prs_pheno,prs_outcome[,c('eid','sex','age')],by.x=c("IID"),by.y='eid') 
  df_prs_pheno<-df_prs_pheno %>% subset(.,outcome<outcome_threshold)
  df_prs_pheno$score_factor<- ifelse(df_prs_pheno$score<quantile(df_prs_pheno$score,0.2
  ),0,ifelse(
    df_prs_pheno$score<quantile(df_prs_pheno$score,0.8
    ),1,2
  )) %>% factor()
  model <- lm(outcome ~ score_factor +sex+age, data = df_prs_pheno )
  return(list(model,df_prs_pheno))
  
}







sp_score<-'./sum_score_SPd18.txt'



summary(get_model(sp_score,11.1)[[1]])

# ==================== 箱线图（含统计比较与数值标注） ====================
# 首先计算每组的均值、四分位数等，用于标注
prs_outcome<- get_model(sp_score,11.1)[[2]]
summary_stats <- prs_outcome %>%
  group_by(score_factor) %>%
  summarise(
    N = n(),
    Mean = mean(outcome, na.rm = TRUE),
    SD = sd(outcome, na.rm = TRUE),
    Q1 = quantile(outcome, 0.25, na.rm = TRUE),
    Median = median(outcome, na.rm = TRUE),
    Max=max(outcome,na.rm=TRUE),
    Min=min(outcome,na.rm=TRUE),
    Q3 = quantile(outcome, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))  # 保留两位小数

# 生成标注文本，例如 "Mean = 12.3\nQ1 = 10.1\nMedian = 11.5\nQ3 = 13.8"
annotation_labels <- summary_stats %>%
  mutate(
    label = paste0(
      "Mean = ", Mean, "\n",
      "Q1 = ", Q1, "\n",
      "Median = ", Median, "\n",
      "Q3 = ", Q3
    )
  ) %>%
  select(score_factor, label)

# 确定标注的 y 坐标（放在每个箱子的上方，略高于最大值）
max_y <- max(prs_outcome$outcome, na.rm = TRUE)
y_offset <- 0.1 * (max_y - min(prs_outcome$outcome, na.rm = TRUE))
y_positions <- summary_stats %>%
  mutate(y_pos = Q3 + y_offset)  # 放在第三四分位数上方

# 合并标注信息
annotations <- left_join(y_positions, annotation_labels, by = "score_factor")

# 箱线图绘制
p <- ggplot(prs_outcome, aes(x = score_factor, y = outcome, fill = score_factor)) +
  geom_violin(alpha = 0.5, outlier.shape = NA, width = 0.5) +   # 箱线图，不显示离群点（避免与抖动点重叠）
  geom_boxplot(width=0.1,outliers = FALSE, box.color = 'black', fill='white',lwd=0.5)+
  #geom_jitter(width = 0.15, alpha = 0.3, size = 1.5, color = "gray30") +  # 原始数据点
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "darkgrey") +  # 均值点
  scale_fill_viridis(discrete = TRUE, begin = 0.2, end = 0.8, name = "Score") +        # 科研配色
  labs(
    title = "Group comparisons",
    x = "Polygenic score group",
    y = "Fasting glucose"
  ) +  
  scale_x_discrete(
    labels = c("Low","Medium","High")
  )+
  theme_classic(base_size = 10) +                     # 简洁主题
  theme(
    legend.position = "none",                          # 删除图例（颜色与x轴分组一致）
    plot.title = element_text(hjust =0.7, face = "bold",size=8), # 主标题
    axis.title = element_text(size = 8), # 坐标轴标题
    axis.text = element_text(size = 8) # 坐标轴刻度
    
  )

# 添加两两比较的显著性标记
# 获取所有可能的组间比较
score_levels <- levels(prs_outcome$score_factor)
comparisons_list <- combn(score_levels, 2, simplify = FALSE)

# 使用 ggpubr 的 stat_compare_means 添加 p 值（wilcox.test 或 t.test）
p <- p + stat_compare_means(
  comparisons = comparisons_list,
  method = "wilcox.test",            # 可选 t.test 或 wilcox.test
  label = "p.signif",                 # 显示显著性符号（*，**，***，ns）
  symnum.args = list(
    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    symbols = c("***", "**", "*", "ns")
  ),
  tip.length = 0.02,
  bracket.size = 0.3,
  step.increase = 0.1,
  label.y =10.8
)

# 添加统计量标注（均值、四分位数）
#annotations$score=levels(annotations$score) %>% as.numeric()
p <- p + geom_text(
  data = annotations,
  aes(x = score_factor, y = y_pos+3, label = label),
  size = 2.9,
  hjust = -0.1,
  color = "black"
)
p
plot_list[['box_sp_11']]<-p

library(plotRCS)


# 线性回归模型的 RCS 曲线
rcsplot(data = prs_outcome %>% subset(.,outcome<7),
        outcome = "outcome",
        exposure = "score",
        covariates = c( "age",'sex'))
plot(prs_outcome$score,prs_outcome$outcome)


library(quantreg)
# 拟合多个分位数
taus <- seq(0.1,0.9,0.1)
qr_all <- rq(outcome ~ score + sex+age,data = prs_outcome %>% subset(outcome<11.101), tau = taus)
qr_sum <- summary(qr_all, se = "boot") 

# 提取 score 的信息
beta_score <- sapply(qr_sum, function(x) x$coefficients["score", 1])  %>% as.numeric() # 系数
se <- sapply(qr_sum, function(x) x$coefficients["score", 2])  %>% as.numeric() # se
upper_score <- beta_score +1.96*se # 上界
lower_score<-  beta_score -1.96*se # 上界
pval_score  <- sapply(qr_sum, function(x) x$coefficients["score", 4])  %>% as.numeric()# P值

# 创建表格
coef_score <- data.frame(
  tau = taus,
  beta = beta_score,
  lower = lower_score,
  upper = upper_score,
  p.value = pval_score
)

# 可选：添加显著性标记
coef_score$sig <- ifelse(coef_score$p.value < 0.05, "*", "")
coef_score$sig2 <- ifelse(coef_score$p.value < 0.01, "**", 
                          ifelse(coef_score$p.value < 0.05, "*", ""))

# 查看表格
print(coef_score)
library(ggplot2)

# 添加一个格式化 P 值的列
coef_score$p_label <- paste0("P = ", sprintf('%.2e',coef_score$p.value) )

ggplot(coef_score, aes(x = tau, y = beta)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(size = 2) +
  geom_text(aes(label = p_label, y = beta + 0.05 * max(beta)),  # 向上偏移一点
            size = 3, hjust = 0.5, vjust = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Quantile (τ)", y = "Effect of score") +
  theme_bw()


# 假设 coef_score 已存在，包含 tau, beta, lower, upper, p.value

# 创建两个标签列
coef_score$beta_label <- sprintf("%.2f [%.2f, %.2f]", 
                                 coef_score$beta, 
                                 coef_score$lower, 
                                 coef_score$upper)
coef_score$p_label <- sprintf("P = %.2e", coef_score$p.value)
coef_score$tau_label <- paste0(round(coef_score$tau * 100), "th")
# 计算垂直偏移量（基于 beta 的范围）
y_range <- diff(range(coef_score$beta, na.rm = TRUE))
offset_up <- 0.1 * y_range   # 向上偏移量
offset_down <- 0.1 * y_range # 向下偏移量，如果怕重叠可稍微减小

quant_plot_sp<-ggplot(coef_score, aes(x = tau, y = beta)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(size = 2) +
  
  # 上标签：beta [lower, upper]
  geom_text(aes(label = beta_label, y = beta + offset_up), 
            size = 2.1, hjust = 0.5, vjust = 0) +
  # 下标签：P 值
  geom_text(aes(label = p_label, y = beta - offset_down+0.01), 
            size = 2.3, hjust = 0.5, vjust = 0.3) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = coef_score$tau, 
                     labels = coef_score$tau_label) +
  labs(x = "Quantile of fasting glucose (τ)", y = "Effect (beta)") +
  theme_bw()
quant_plot_sp
plot_list[['quant_sp_11']]<-quant_plot_sp





# 拟合模型
lm_fit <- lm(outcome ~ score_factor + sex+age, data = prs_outcome)
# 提取系数和置信区间
sort_lmdata<- function(lm_fit){
  library(broom)
  tidy_lm <- tidy(lm_fit, conf.int = TRUE, conf.level = 0.95)
  
  # 提取 score_factor 的系数（自动排除参考组 low）
  tidy_score <- tidy_lm[grep("score_factor", tidy_lm$term), ]
  # 修正 term 名称：去掉 "score_factor" 前缀
  tidy_score$term <- gsub("score_factor", "", tidy_score$term)
  
  # 创建参考组空行
  ref_row <- data.frame(
    term = "low (ref)",
    estimate = NA,
    conf.low = NA,
    conf.high = NA,
    p.value = NA
  )
  
  # 合并，确保参考组在第一行（后面会反转顺序）
  tidy_score$term<-c('medium','high')
  forest_data <- rbind(ref_row, tidy_score[, c("term", "estimate", "conf.low", "conf.high", "p.value")])
  
  # 设置因子顺序：从上到下为 high, medium, low (ref)
  forest_data$term <- factor(forest_data$term, 
                             levels = c("high", "medium", "low (ref)"))
  
  # 为每个非参考组生成文本标签：beta (95% CI)
  forest_data$beta_ci <- ifelse(!is.na(forest_data$estimate),
                                sprintf("%.2f (%.2f, %.2f)", 
                                        forest_data$estimate, 
                                        forest_data$conf.low, 
                                        forest_data$conf.high),
                                "")
  
  forest_data$p_label <- ifelse(!is.na(forest_data$p.value),
                                sprintf("P = %.2e", forest_data$p.value),
                                "")
  return (forest_data)
}
whole_lm=sort_lmdata(lm_fit)
prediabete <- lm(outcome ~ score_factor + sex+age, data = prs_outcome %>% subset(outcome<7))
pre_lm=sort_lmdata(prediabete)
health_lm<- lm(outcome ~ score_factor + sex+age, data = prs_outcome %>% subset(outcome< 6.1)) %>% sort_lmdata()

whole_lm$group<- 'fasting glucose < 11.1'
pre_lm$group <- 'fasting glucose < 7'
health_lm$group<- 'fasting glucose < 6.1'
forest_data=rbind(whole_lm,pre_lm,health_lm)

# 确定 x 轴的最大值（用于定位右侧文本）
x_max <- forest_data$conf.high #max(forest_data$conf.high, na.rm = TRUE)
x_range <- diff(range(c(0, x_max), na.rm = TRUE))
# 右侧文本的 x 坐标：x_max + 0.2 * x_range
x_right <- x_max + 0.05 * x_range


forest_data$group<- factor(forest_data$group,levels=c('fasting glucose < 11.1','fasting glucose < 7','fasting glucose < 6.1'))



# 绘制森林图，并添加右侧文本
popFg_group_forest<-ggplot(forest_data, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = p.value < 0.05), size = 1.5, na.rm = TRUE) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = p.value < 0.05),
                 height = 0.2, na.rm = TRUE) +
  # 右侧文本：β (95% CI)
  geom_text(aes(x = x_right, label = beta_ci), 
            hjust = 0-0.2, vjust=0,size = 2.5, na.rm = TRUE) +
  # 右侧文本：P值（放在 CI 文本右边一点）
  geom_text(aes(x = x_right + 0.2 * x_range, label = p_label), 
            hjust=0-0.4,
            vjust =2-0.2, size = 2.5, na.rm = TRUE) +
  scale_color_manual(values = c("black", "red"), 
                     name = "P-value",
                     labels = c("P ≥ 0.05", "P < 0.05"),
                     na.translate = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.6))) +  # 右侧留更多空间
  labs(x = "Effect (beta,95% CI)", y = "PGS level") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 11))+
  facet_wrap(vars(group), nrow=3)
popFg_group_forest
plot_list[['popForest']]<-popFg_group_forest



wc_score='./sum_score_WC.txt'
prs_outcome<- get_model(wc_score,11.101)[[2]]
taus <- seq(0.1,0.9,0.1)
qr_all <- rq(outcome ~ score + sex+age,data = prs_outcome %>% subset(outcome<11.101), tau = taus)
qr_sum <- summary(qr_all, se = "boot") 

# 提取 score 的信息
beta_score <- sapply(qr_sum, function(x) x$coefficients["score", 1])  %>% as.numeric() # 系数
se <- sapply(qr_sum, function(x) x$coefficients["score", 2])  %>% as.numeric() # se
upper_score <- beta_score +1.96*se # 上界
lower_score<-  beta_score -1.96*se # 上界
pval_score  <- sapply(qr_sum, function(x) x$coefficients["score", 4])  %>% as.numeric()# P值

# 创建表格
coef_score <- data.frame(
  tau = taus,
  beta = beta_score,
  lower = lower_score,
  upper = upper_score,
  p.value = pval_score
)

# 可选：添加显著性标记
coef_score$sig <- ifelse(coef_score$p.value < 0.05, "*", "")
coef_score$sig2 <- ifelse(coef_score$p.value < 0.01, "**", 
                          ifelse(coef_score$p.value < 0.05, "*", ""))

# 查看表格
print(coef_score)
library(ggplot2)

# 添加一个格式化 P 值的列
coef_score$p_label <- paste0("P = ", sprintf('%.2e',coef_score$p.value) )

ggplot(coef_score, aes(x = tau, y = beta)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(size = 2) +
  geom_text(aes(label = p_label, y = beta + 0.05 * max(beta)),  # 向上偏移一点
            size = 3, hjust = 0.5, vjust = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Quantile (τ)", y = "Effect of score") +
  theme_bw()


# 假设 coef_score 已存在，包含 tau, beta, lower, upper, p.value

# 创建两个标签列
coef_score$beta_label <- sprintf("%.2f [%.2f, %.2f]", 
                                 coef_score$beta, 
                                 coef_score$lower, 
                                 coef_score$upper)
coef_score$p_label <- sprintf("P = %.2e", coef_score$p.value)
coef_score$tau_label <- paste0(round(coef_score$tau * 100), "th")
# 计算垂直偏移量（基于 beta 的范围）
y_range <- diff(range(coef_score$beta, na.rm = TRUE))
offset_up <- 0.1 * y_range   # 向上偏移量
offset_down <- 0.1 * y_range # 向下偏移量，如果怕重叠可稍微减小

ggplot(coef_score, aes(x = tau, y = beta)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(size = 2) +
  
  # 上标签：beta [lower, upper]
  geom_text(aes(label = beta_label, y = beta + offset_up), 
            size = 2.1, hjust = 0.5, vjust = 0) +
  # 下标签：P 值
  geom_text(aes(label = p_label, y = beta - offset_down-0.2), 
            size = 2.3, hjust = 0.5, vjust = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = coef_score$tau, 
                     labels = coef_score$tau_label) +
  labs(x = "Quantile of fasting glucose (τ)", y = "Effect (beta)") +
  theme_bw()

rcsplot(data = prs_outcome %>% subset(.,outcome<11.1),
        outcome = "outcome",
        exposure = "score",
        covariates = c( "age",'sex'))

wc_score='./sum_score_WC.txt'
prs_outcome<- get_model(wc_score,11.101)[[2]]
taus <- seq(0.1,0.9,0.1)
qr_all <- rq(outcome ~ score + sex+age,data = prs_outcome %>% subset(outcome<11.101), tau = taus)
qr_sum <- summary(qr_all, se = "boot") 

# 提取 score 的信息
beta_score <- sapply(qr_sum, function(x) x$coefficients["score", 1])  %>% as.numeric() # 系数
se <- sapply(qr_sum, function(x) x$coefficients["score", 2])  %>% as.numeric() # se
upper_score <- beta_score +1.96*se # 上界
lower_score<-  beta_score -1.96*se # 上界
pval_score  <- sapply(qr_sum, function(x) x$coefficients["score", 4])  %>% as.numeric()# P值

# 创建表格
coef_score <- data.frame(
  tau = taus,
  beta = beta_score,
  lower = lower_score,
  upper = upper_score,
  p.value = pval_score
)

# 可选：添加显著性标记
coef_score$sig <- ifelse(coef_score$p.value < 0.05, "*", "")
coef_score$sig2 <- ifelse(coef_score$p.value < 0.01, "**", 
                          ifelse(coef_score$p.value < 0.05, "*", ""))

# 查看表格
print(coef_score)
library(ggplot2)

# 添加一个格式化 P 值的列
coef_score$p_label <- paste0("P = ", sprintf('%.2e',coef_score$p.value) )

ggplot(coef_score, aes(x = tau, y = beta)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(size = 2) +
  geom_text(aes(label = p_label, y = beta + 0.05 * max(beta)),  # 向上偏移一点
            size = 3, hjust = 0.5, vjust = 0) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Quantile (τ)", y = "Effect of score") +
  theme_bw()


# 假设 coef_score 已存在，包含 tau, beta, lower, upper, p.value

# 创建两个标签列
coef_score$beta_label <- sprintf("%.2f [%.2f, %.2f]", 
                                 coef_score$beta, 
                                 coef_score$lower, 
                                 coef_score$upper)
coef_score$p_label <- sprintf("P = %.2e", coef_score$p.value)
coef_score$tau_label <- paste0(round(coef_score$tau * 100), "th")
# 计算垂直偏移量（基于 beta 的范围）
y_range <- diff(range(coef_score$beta, na.rm = TRUE))
offset_up <- 0.1 * y_range   # 向上偏移量
offset_down <- 0.1 * y_range # 向下偏移量，如果怕重叠可稍微减小

quant_wc5e6<- ggplot(coef_score, aes(x = tau, y = beta)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point(size = 2) +
  
  # 上标签：beta [lower, upper]
  geom_text(aes(label = beta_label, y = beta + offset_up+0.1), 
            size = 2.1, hjust = 0.5, vjust = 0) +
  # 下标签：P 值
  geom_text(aes(label = p_label, y = beta - offset_down-0.07), 
            size = 2.3, hjust = 0.5, vjust = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = coef_score$tau, 
                     labels = coef_score$tau_label) +
  labs(x = "Quantile of fasting glucose (τ)", y = "Effect (beta)") +
  theme_bw()
quant_wc5e6
plot_list[['quant_wc']]<-quant_wc5e6


rcs_wc<- rcsplot(data = prs_outcome %>% subset(.,outcome<11.1),
        outcome = "outcome",
        exposure = "score",
        covariates = c( "age",'sex'))+
  theme(
    axis.title.x=element_text(family = "sans"),
    axis.title.y=element_text(family = "sans"),
    axis.text.x = element_text(family = "sans"),  # x 轴刻度文字
    axis.text.y = element_text(family = "sans")   # y 轴刻度文字
  )
rcs_wc
plot_list[['rcs_wc']]<-rcs_wc

library(cowplot)
p1=plot_grid(plot_list[[1]],plot_list[[3]],ncol=2,rel_widths=c(1.5, 1),labels=c('a','b'))
p2=plot_grid(plot_list[[2]],plot_list[[5]],plot_list[[4]],nrow=3,labels=c('c','d','e'))
p_all=plot_grid(p1,p2,nrow=2,rel_heights =c(1.2,2))
ggsave('pAll.png',p_all,height=12,width=8)
plot_grid(plotlist = plot_list,ncol=2,nrow=3,rel_widths=c(1, 2),rel_heights = c(1, 2),labels=c('a','b'))
