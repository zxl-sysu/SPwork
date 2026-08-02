import pandas as pd
import os
import sys
from tqdm import tqdm as tqdm
import multiprocessing


def pre_qc(raw_gwas):
    pass

    return raw_gwas
def match_a_chr_pos(match,chr):
    snp_remove=[]
    flip=[]
    swap=[]
    flip_swap=[]
    if not match.empty:
        #match=pd.merge(match,frq.loc[:,['rsID','A1_frq','A2_frq','MAF']])
        print('translate complete in chr%s'%chr)
        
        complete_dict={'A':'T','G':'C','T':'A','C':'G'}
        match['C.ea']=match['effect_allele'].apply(lambda  x: complete_dict[x])
        match['C.oa']=match['other_allele'].apply(lambda  x: complete_dict[x])


        print('complete done, matching swap flip in chr%s'%chr)
        match['swap']=(match['effect_allele']==match['A2'] ) & (match['other_allele']==match['A1'] ) 
        match['flip']= (match['C.ea']==match['A1'] ) & (match['C.oa']==match['A2'] ) 
        match['flip_swap']= (match['C.ea']==match['A2'] ) & (match['C.oa']==match['A1'] ) 
        

        com=[{'A','T'},{'G','C'}]
        match['complete']=match.apply(lambda x: set(x.loc[
            ['effect_allele','other_allele']].to_list()
            ) in  com, axis=1)
        
        


        match['pg_gene']=match.apply(lambda x: set(x.loc[
            ['effect_allele','other_allele']].to_list()
            ),  axis=1)
        match['bim_gene']=match.apply(lambda x: set(x.loc[
            ['A1','A2']].to_list()
            ),  axis=1)
        match['gene_match']=match['bim_gene']==match['pg_gene']

        print('match done, saving complete swap flip in chr%s'%chr)

        complete=match[match.complete]
        
        mismatch=match[(~match['gene_match']) & (~match['flip']) & (~match['flip_swap']) &(~ match['complete'])]
        if not mismatch.empty:
            #display('mis_gene_match',mismatch)
            
            snp_remove.extend(mismatch.rsID.to_list())
        if not complete.empty:
            #display('complete',complete)
            snp_remove.extend(complete.rsID.to_list())
        
        #display(all(match.A1==match.A1_frq))

        if not match[match['flip']].empty:
            flip+=match[match['flip']].rsID.to_list()
        if not match[match['swap']].empty:
            swap+=match[match['swap']].rsID.to_list()
            #swap=match[match['swap']]
            #display( (swap['EAF']+swap['MAF']).hist() )
        if not match[match['flip_swap']].empty:
            flip_swap+=match[match['flip_swap']].rsID.to_list()
    match=None
    return([snp_remove,flip,swap,flip_swap])

if __name__=='__main__':

    gwas_dir= 'myquery.txt' ##改成base gwas的路径
    

    save_dir='.'
    os.makedirs(save_dir,exist_ok=True)
    print('reading BASE DATA from ...%s' %gwas_dir)






    pg1_raw=pd.read_csv(gwas_dir ,sep='\t') #若分隔符不是'\t', 改为对应的分隔符
    print('read done!...%s' % gwas_dir)
    pg1_raw=pre_qc(pg1_raw)

    ## rsID改成 snp名所在的列名
    pg1=pg1_raw[~pd.isna(pg1_raw['SNP']) ]
    pg1=pg1[pg1['Gene']=='TMEM258']
    
    ## 左侧是base gwas中相应列的列名，右侧是标准列名，只需改左侧，如rsID代表SNP名，若对应的列是SNP，改为SNP

    rename_dict={
        'SNP':'rsID',
    'Chr':'chr_name',
    'BP':'chr_position',
    'A1':'effect_allele',
    'A2':'other_allele',
    'p':'P-value',
    'Freq':'EAF',
    'b':'effect_weight'}
    
    pg1=pg1.loc[:,rename_dict.keys()].rename(rename_dict,axis=1)
    #pg1['chr_name']=pg1.chr_name.str.replace('chr','').astype('int32')
    pg1['effect_allele']=pg1['effect_allele'].str.upper()
    pg1['other_allele']=pg1['other_allele'].str.upper()
    
    print(pg1.head())

    def tran_chr(x):
        if isinstance(x,str):
            try:
                x2=int(x)
            except:
                x2=23
        else:
            x2=x
        return x2
    #pg1.chr_name=pg1.chr_name.apply(lambda x:  tran_chr(x)   ) 


    ## remove Duplicate SNPs
    pg1=pg1[~pg1.rsID.duplicated()]

    # remove indels
    pg1=pg1[~((pg1.effect_allele.str.len()>1) | (pg1.other_allele.str.len()>1))]


    print('start identify mismathcing filp swap...' )
    # identify mismatching SNPs
    target_file="./impute1420"
    ##这是1420人数据路径，若要使用UKB数据，改成UKB路径

    snp_remove=[]
    flip=[]
    swap=[]
    flip_swap=[]
    to_be_modify={}

    ##return([snp_remove,flip,swap,flip_swap])
    snp_tobe_change=[]
    def cbk(x):
        snp_tobe_change.append(x)

    pool=multiprocessing.Pool()

    for i in tqdm(range(1,23),total=22,desc='chr'):
        clp_score=pg1[pg1.chr_name==i].copy(deep=True)
        bim_file='%s/chr%s.bim'%(target_file,i)
        bim=pd.read_csv(bim_file,sep='\t',header=None)
        bim.columns=['chr_name','rsID','mol','chr_position_bim','A1','A2']

        

        match=pd.merge(clp_score.loc[:,['rsID','effect_allele','other_allele','EAF'] ],
        bim.loc[:,['rsID','chr_position_bim','A1','A2']],how='inner')

        pool.apply_async(match_a_chr_pos,args=(match,i),callback=cbk)
        match=None
    pool.close()
    pool.join()

    print('match done, start to modify')
    for i in snp_tobe_change:
        snp_remove.extend(i[0])
        flip.extend(i[1])
        swap.extend(i[2])
        flip_swap.extend(i[3])
        print('remove %s flip %s swap %s flip_swap %s'%tuple([len(j) for j in i ])  )

    # flip+complete=swap   G C--filp= C G swap= C G 


    to_be_modify['flip']=flip
    to_be_modify['swap']=swap
    to_be_modify['flip_swap']=flip_swap
    complete_dict={'A':'T','G':'C','T':'A','C':'G'}

    # remove mismatch snp and ambigulous snp (complementary ,A-T and G-C)
    pg1=pg1[~pg1.rsID.isin(snp_remove)].copy(deep=True)


    print('%s snp to be flipped' % len( pg1[pg1.rsID.isin(to_be_modify['flip'])] ))
    pg1.loc[pg1.rsID.isin(to_be_modify['flip']),'effect_allele']=pg1[pg1.rsID.isin(to_be_modify['flip'])].effect_allele.apply(lambda x:complete_dict[x])
    pg1.loc[pg1.rsID.isin(to_be_modify['flip']),'other_allele']=pg1[pg1.rsID.isin(to_be_modify['flip'])].other_allele.apply(lambda x:complete_dict[x])

    print('%s snp to be flipped_swapped' % len( pg1[pg1.rsID.isin(to_be_modify['flip_swap'])] ))
    pg1.loc[pg1.rsID.isin(to_be_modify['flip_swap']),'effect_allele']=pg1[pg1.rsID.isin(to_be_modify['flip_swap'])].effect_allele.apply(lambda x:complete_dict[x])
    pg1.loc[pg1.rsID.isin(to_be_modify['flip_swap']),'other_allele']=pg1[pg1.rsID.isin(to_be_modify['flip_swap'])].other_allele.apply(lambda x:complete_dict[x])

    to_be_swap=pg1.loc[pg1.rsID.isin(to_be_modify['flip_swap']) | pg1.rsID.isin(to_be_modify['swap']) ].copy(deep=True)

    print('%s snp to be swapped' % len( to_be_swap))
    pg1.loc[pg1.rsID.isin(to_be_modify['flip_swap'])
            | pg1.rsID.isin(to_be_modify['swap']),'effect_allele' ]=to_be_swap['other_allele']
    pg1.loc[pg1.rsID.isin(to_be_modify['flip_swap'])
            | pg1.rsID.isin(to_be_modify['swap']),'other_allele' ]=to_be_swap['effect_allele']
    pg1.loc[pg1.rsID.isin(to_be_modify['flip_swap'])
            | pg1.rsID.isin(to_be_modify['swap']),'effect_weight' ]=-to_be_swap['effect_weight']
    pg1.loc[pg1.rsID.isin(to_be_modify['flip_swap'])
            | pg1.rsID.isin(to_be_modify['swap']),'EAF' ]=1-to_be_swap['EAF']

    print('identify mismathcing filp swap done!... ' )

    ## Sample overlap
    ## it shoulb be aviod by selecting GWAS and target data from different populations

    # remove snps with MAF<.0.01
    pg1=pg1[(pg1.EAF>=0.01) & (pg1.EAF<=0.99)]


    print('saving qc data to qc_base.txt......' )
    pg1_dir='%s/qc_base.txt' % save_dir   
    pg1.to_csv(  pg1_dir ,sep=' ',index=False)

    ##write filped records


    print('saving flipped swapped and mismatch snp (have been removed from qc data) to flip_swap_mismatch.xlsx ......' )
    change_dir='%s/flip_swap_mismatch.xlsx' % save_dir
    with pd.ExcelWriter(change_dir) as f:
        for i,j in to_be_modify.items():
            pd.DataFrame(j,columns=[i]).to_excel(f,sheet_name=i,index=False)
        pd.DataFrame(snp_remove,columns=['mismatch']).to_excel(f,sheet_name='mismatch',index=False)

    print('==============DONE QC of base data=========')





