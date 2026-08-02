import os
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats
from multiprocessing import Pool
import sys


from tqdm import tqdm as tqdm

def func_cmd(cmd):
    os.system(cmd)
    #print(cmd)
if __name__=='__main__':

    ##这是200人数据路径，若要使用UKB数据，改成UKB路径
    target_file=r'./impute1420'

    PRS_dir='.'
    disease_name='SPd18' ##改成表型的名称，用于命名结果文件
    os.chdir(PRS_dir)
    print('PRS dir = %s' % PRS_dir)
    print('changing working dir to %s'% os.getcwd() )

    pt=[5e-8,5e-7,5e-6,1e-5] #若有希望运行的其它P阈值，可修改 5e-1,5e-2,5e-3,5e-4,5e-5,5e-6,5e-7,5e-8
    pt.sort()
    print(pt)
    clump_kb=1000
    clump_r2_list=[0.1,0.01,0.001,0.99] #若希望运行其它clump r2阈值，可修改
    print(clump_r2_list)
    clump_dirs=['clump_1000_r2_%s'%i for i in clump_r2_list]
    for clp_dir in clump_dirs:
        if os.path.exists(clp_dir):
            #os.system('rm -r %s' % clp_dir)
            pass

    for clump_r2 in clump_r2_list:

        #clump_r2=0.001

        trait='%s' % (disease_name)

        clp_dir_name='clump_%s_r2_%s' % (clump_kb,clump_r2)

        score=pd.read_csv(r'./qc_base.txt',sep=r'\s+')


        pool_clp=Pool(len(pt))
        for p in tqdm(pt,desc='clumping...',total=len(pt)):
            # clumping 
            save_dir='%s/p_%s'%(clp_dir_name,p)
            os.makedirs(save_dir,exist_ok=True)

            old_score_files=[ i for i in os.listdir(save_dir) if i.startswith('%s_prs_chr' % trait) ]
            for i in old_score_files:
                os.remove(os.path.join(save_dir,i))

            ##使用欧洲人群数据时将bfile 文件路径最后一级改成EUR
            clump='plink \
                --bfile  ./1kg.v3/EAS  \
                --clump-p1 %s \
                --clump-r2 %s \
                --clump-kb %s \
                --clump ./qc_base.txt \
                --clump-snp-field rsID \
                --clump-field P-value \
                --out %s/%s_clumped' % (p,clump_r2,clump_kb,save_dir,trait)
            
            pool_clp.apply_async(func=func_cmd,args=(clump,))

            #os.system(clump)
        pool_clp.close()
        pool_clp.join()

        print('============clump',p,'DNOE============')

        for p in pt:
            save_dir='%s/p_%s'%(clp_dir_name,p)
            print('PRS',save_dir)

            clump_path='%s/%s_clumped.clumped'% (save_dir,trait)

            if os.path.exists(clump_path):
                pass
            else:
                continue
            cho_clp_snp=pd.read_csv(clump_path,sep=r'\s+')

            cho_clp_snp.loc[:,'SNP'].to_csv('%s/%s_snp.txt' %(save_dir,trait) ,index=None,header=None)

            score_clp=score[score.rsID.isin(cho_clp_snp.loc[:,'SNP'].to_list())]


            score_clp.to_csv('%s/%s_score_clp.txt'% (save_dir,trait) ,index=None,sep=' ')

            #1st column is the SNP ID; 4th column is the effective allele information; the 8th column is the effect size estimate;
            score_file='%s/%s_score_clp.txt'% (save_dir,trait)

            pool_prs=Pool(2)
            for i in tqdm(range(1,23),desc='%s_PRS'%p,total=22):
                prs='./plink2.exe \
                    --bfile %s \
                    --chr %s \
                    --score %s 1 4 8 header list-variants cols=+scoresums \
                    --extract %s \
                    --out %s' % ('%s/all_without_miss'%(target_file),i,score_file, 
                    '%s/%s' % (save_dir,'%s_snp.txt' % trait),
                    '%s/%s_chr%s'%(save_dir,'%s_prs' % trait,i))
                os.system(prs)

                #pool_prs.apply_async(func=func_cmd,args=(prs,))
            pool_prs.close()
            pool_prs.join()
            print('============PRS',p,'DNOE============')

            score_dfs=[]
            score_files=[i for i in os.listdir(save_dir) if i.endswith('sscore')]
            if score_files:
                for file in score_files:

                    i=file.split('chr')[-1].split('.sscore')[0]
                    prs_result=pd.read_csv('%s/%s_prs_chr%s.sscore'%(save_dir,trait,i) ,sep=r'\s+')
                    prs_result.columns=['FID'] + prs_result.columns.tolist()[1:]
                    score_df=prs_result.iloc[:,[1,-1]]
                    score_df=score_df.rename({'SCORE1_SUM':'score_chr_%s'%i},axis=1)
                    score_df.set_index('IID',inplace=True)
                    #display(score_df.head())
                    score_dfs.append(score_df)

                all_chr_score=pd.concat(score_dfs,axis=1)

                all_chr_score.to_csv('%s/all_scores.csv'%save_dir,index=False)

                sum_score=all_chr_score.sum(axis=1)
                sum_score.name='score'

                #sum_score.to_csv('%s/sum_score_cho.csv' % save_dir)
                sum_score.to_csv('%s/sum_score_%s.txt' % (save_dir,trait) ,sep='\t')

                print('%s/sum_score_%s.txt')
                ##可以改变直方图bins的数目
                fig=plt.figure()
                plt.hist(sum_score.to_list(),bins=10,range=(sum_score.min(),sum_score.max()))
                plt.savefig('%s/hst50.png'% save_dir )
                plt.close()

                
                

                ks_result=stats.kstest(sum_score.to_list(),'norm',(sum_score.mean(),sum_score.std()))
                ks2=pd.Series(ks_result,index=['statistic','P'])
                ks2.to_csv('%s/KS_test_normal.txt'%save_dir,sep='\t',header=None)

        print('==================Calculating PRS DONE ALL=====================')




            



