

library(data.table)
library(TwoSampleMR)
library(plinkbinr)
library(ieugwasr)
library(phenoscanner)
library(MRPRESSO)
library(openxlsx)
library(stringr)
library(ggplot2)

#gc()

args<-  commandArgs(trailingOnly = T)

save_dir<- args[1]

print(sprintf('reading... %s',args[2]))

target= args[2]   

if_lifted=args[3]

exposure_name<-args[4]
outcome_name<- args[5]

ieu_id<- args[6]

target_gwas<- fread(target) %>% as.data.frame()



print(sprintf('reading done! %s',args[2]))


base_gwas<- target_gwas

sig_p<-5E-8
base_gwas<- base_gwas %>% subset(.,pval.exposure<sig_p)

base_gwas$pheno<- exposure_name
preffix= outcome_name  
#print(base_gwas$pval.outcome< base_gwas$pval.exposure)
min_p_outcome<-min(base_gwas$pval.outcome)


#base_gwas


# beta = log(OR)

exposure<- format_data(
  base_gwas,
  type="exposure",

  snp_col = "SNP.exposure",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
  chr_col = "chr.exposure",
  pos_col = "pos.exposure",
  phenotype_col = 'pheno',
  samplesize_col = 'N.exposure'
  
)



ex_clp_snps<-ld_clump( dplyr::tibble(
  rsid=exposure$SNP, pval=exposure$pval.exposure), 
  plink_bin = get_plink_exe(), bfile = "./1kg.v3/EAS" ,
  clump_kb = 1000,
  clump_r2=0.1)
bool_snp<- exposure$SNP %in% ex_clp_snps$rsid
ex_clp<- exposure[bool_snp,]

out_gwas<- target_gwas
out_gwas$pheno<- outcome_name

outcome<- format_data(
  out_gwas,
  type="outcome",
  snps=ex_clp$SNP,
  snp_col = "SNP.exposure",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "eaf.outcome",
  pval_col = "pval.outcome",
  chr_col = "chr.outcome",
  pos_col = "pos.outcome",
  phenotype_col = 'pheno',
  samplesize_col = 'N.outcome'
  
)

print(nrow(out_gwas))
print(nrow(ex_clp))
print(nrow(outcome))
preffix= outcome_name  
dat_file=sprintf('%s/%s_%s',save_dir,preffix,'mr_datex.csv')
print(dat_file)
write.csv(ex_clp, dat_file, row.names = F)
dat_file=sprintf('%s/%s_%s',save_dir,preffix,'mr_datout.csv')
print(dat_file)
write.csv(outcome, dat_file, row.names = F)

mr_dat<- harmonise_data(ex_clp,outcome)

dat_file=sprintf('%s/%s_%s',save_dir,preffix,'mr_dat.csv')
print(dat_file)
write.csv(mr_dat, dat_file, row.names = F)


out <- directionality_test(mr_dat)
print(out)
dat_file=sprintf('%s/%s_%s',save_dir,preffix,'mr_steiger.csv')
print(dat_file)
write.csv(out, dat_file, row.names = F)
#steiger_sl<-steiger_filtering(mr_dat)


st_filter=steiger_filtering(mr_dat)
mr_dat<-st_filter[st_filter$rsq.exposure>st_filter$rsq.outcome,]

dat_file=sprintf('%s/%s_%s',save_dir,preffix,'mr_dat.csv')
print(dat_file)
write.csv(mr_dat, dat_file, row.names = F)


info_file=sprintf('%s/%s_%s',save_dir,preffix, 'info.txt')
info=sprintf(
    'before clump N= %s \n after clump N= %s \n outcome N= %s \n harmonise N= %s',
    dim(exposure)[1],
    dim(ex_clp)[1],
    dim(outcome)[1],
    dim(mr_dat)[1])
write.csv(
  info,  
  info_file

)


mr(mr_dat)








do_mr_3method2<-function(ex_outcome){
  MR_result_list=list()
  MR_result_list[['raw_ex']]<-exposure
  MR_result_list[['SNP']]<-ex_outcome
  
  mr_results<-generate_odds_ratios(mr(ex_outcome,method_list = 
    subset(mr_method_list(), use_by_default)$obj[1:4] ) ) 
  print(mr_results)
  MR_result_list[['MR' ]]<-mr_results
  
  MR_result_list[['hete']]<- mr_heterogeneity(ex_outcome,method_list ='mr_egger_regression' )
  
  MR_result_list[['pleio']]<- mr_pleiotropy_test(ex_outcome)
  return(MR_result_list)
}

mr_res<- do_mr_3method2(mr_dat)
mr_res[['steiger']]<- out
mr_res[['steiger_filter']]<-st_filter

#print(sprintf('MRpresso... %s',args[2]))
try(
  { raise
    presso_dat<- as.data.frame(mr_dat)  %>% subset(.,mr_keep)
    presso<-  mr_presso(
    BetaOutcome = 'beta.outcome',
    BetaExposure = 'beta.exposure',
    SdOutcome = 'se.outcome',
    SdExposure ='se.exposure',
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    data = presso_dat ,
    NbDistribution = nrow(presso_dat)/0.05 ,
    SignifThreshold = 0.05) ##[[2]][[1]]$Pvalue
  
  
  mr_res[['presso']]<- presso$`Main MR results`
  print(mr_res[['presso']])
  }
)
#print(sprintf('MRpresso done! %s',args[2]))

add_format<- function(mr_plot_obj) {
  p1<- mr_plot_obj [[1]]
  p1<- p1 + theme(text = element_text(family = "Times New Roman"))
  return(p1)
  }

res_single <- mr_singlesnp(mr_dat)
mr_res[['leave_one_out']]<- mr_leaveoneout(mr_dat)
mr_res[['single']]<-res_single
res_name<- sprintf('%s/%s_%s',save_dir,preffix,'MR_res.xlsx')
write.xlsx(mr_res,res_name)

stop()
lve_plot<- mr_leaveoneout_plot(  mr_res[['leave_one_out']] )  
plot_name<- sprintf('%s/%s_%s',save_dir,preffix,'leave_one_out.png')
row_n=dim(mr_res[['leave_one_out']])[1]

ggsave(plot_name,add_format(lve_plot),height = row_n*5,limitsize = F,units = 'mm')


sct <- mr_scatter_plot (mr_res[['MR']], mr_dat)
plot_name<- sprintf('%s/%s_%s',save_dir,preffix,'scatter.png')
ggsave(plot_name,add_format(sct),limitsize = T)


forest_single<- mr_forest_plot(res_single)
row_n=dim(res_single)[1]
plot_name<- sprintf('%s/%s_%s',save_dir,preffix,'forest_single.png')
ggsave(plot_name,add_format(forest_single),height = row_n*5,limitsize = F,units = 'mm')

funnel<- mr_funnel_plot(res_single)
plot_name<- sprintf('%s/%s_%s',save_dir,preffix,'funnel.png' )
ggsave(plot_name,add_format(funnel),limitsize = F)



gc()
