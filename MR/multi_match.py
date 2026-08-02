import pandas as pd
import os
import multiprocessing as mtp

def get_a_ex_a_out(ex_name,ex_df,out_name,outdf,save_dir):

    print(ex_name,'doing outcome %s ...' % out_name)
    
    
    print('doing exposure %s ...' % ex_name)
    #print(ex_df.head() )
    
    df_out=outdf

    
    
    
    
    
    

    match_ex_out=pd.merge(ex_df,df_out,left_on=[
        'chr.exposure','pos.exposure'],right_on=[
            'chr.outcome','pos.outcome'],how='inner')
    print(match_ex_out)
    if not match_ex_out.empty:
        save_path='%s/%s/%s' % (save_dir,ex_name,out_name)
        os.makedirs(save_path,exist_ok=True)

        save_file='%s/ex_out_sig.txt' % (save_path)

        match_ex_out.to_csv(save_file,sep='\t',index=None)

        #print('exposure %s outcome %s done!' % (ex_name,out_name))
        return('exposure %s outcome %s done!' % (ex_name,out_name))
    else:
        #print('exposure %s outcome %s no snp matched' % (ex_name,out_name))
        return('exposure %s outcome %s no snp matched' % (ex_name,out_name)) 


if __name__=='__main__':
    

    save_dir='match_ex_out'

    # before: we extract the exposure data with P < 1E-5 to this path
    ex_dir_all='./Sphingolipids/sig'

    ex_list=os.listdir(ex_dir_all)
    ex_dict={}
    for ex_name_raw in ex_list:
        ex_name=ex_name_raw.split('.txt')[0]
        ex_path='%s/%s' % (ex_dir_all,ex_name_raw)
        ex_df=pd.read_csv(ex_path,sep='\t')

        ex_df=ex_df.rename(
            {'CHR':'chr.exposure',
            'POS':'pos.exposure',
            'EFF_ALLELE':'effect_allele.exposure',
            'NONEFF_ALLELE':'other_allele.exposure',
            'N':'N.exposure',
            'EFF_ALLELE_FREQ':'eaf.exposure',
            'MARKERNAME':'SNP.exposure',
            'BETA':'beta.exposure',
            'SEBETA':'se.exposure',
            'P':'pval.exposure'
            },axis=1
        )
        ex_df['effect_allele.exposure']=ex_df['effect_allele.exposure'].str.upper()
        ex_df['other_allele.exposure']=ex_df['other_allele.exposure'].str.upper()
        ex_df=ex_df[(ex_df['effect_allele.exposure'].str.len()==1) & (
            ex_df['other_allele.exposure'].str.len()==1)    ]
        
        ex_df=ex_df[(ex_df['eaf.exposure']<0.99) &( ex_df['eaf.exposure']> 0.01)]
        
        ex_df=ex_df.drop_duplicates(subset=['chr.exposure','pos.exposure']).copy(deep=True)
        ex_dict[ex_name]=ex_df
    ex_snps=[list(zip(j['chr.exposure'].to_list(),j['pos.exposure'].to_list())) for i,j in ex_dict.items()]
    #print([j for i in ex_snps for j in i])
    ex_snps=['%s:%s'%(j[0],j[1]) for i in ex_snps for j in i]



    out_dict={}
    out_dict['T2D']="./Suzuki.Nature2024.T2DGGI.EAS.sumstats.zip"
    out_df_dict={}
    for out_name,out_path in out_dict.items():
        outdf=pd.read_csv(out_path,sep='\t')
        
        
        outdf=outdf.rename(
            {#'ID':'SNP.outcome',
            'Chromsome':'chr.outcome',
            'Position':'pos.outcome',
            'EffectAllele':'effect_allele.outcome',
            'NonEffectAllele':'other_allele.outcome',
            
            'EAF':'eaf.outcome',
            'Beta':'beta.outcome',
            'SE':'se.outcome',
            'Pval':'pval.outcome'
            },axis=1
        )
        outdf['N.outcome']= outdf['Ncases']+outdf['Ncontrols']	
        print(outdf.head())
        outdf=outdf.drop_duplicates(subset=['chr.outcome','pos.outcome']).copy(deep=True)

        outdf['cptid.outcome']=outdf['chr.outcome'].astype('str') + ':' + outdf['pos.outcome'].astype('str')
            
            
        print(ex_snps[:5])
        outdf=outdf[outdf['cptid.outcome'].isin(ex_snps)]
        out_df_dict[out_name]=outdf
    pool=mtp.Pool()
    ex_dir_all='/data/user1/zixin_work/QTL_SMR/agen_metabolite/Sphingolipids/sig'

    def cbk(x):
        print(x)

    for ex_name, ex_df in ex_dict.items():
        for out_name,outdf in out_df_dict.items():         
            #print('doing outcome %s ...' % out_name)
            ex_path='%s/%s' % (ex_dir_all,ex_name_raw)
            pool.apply_async(get_a_ex_a_out,args=(ex_name,ex_df,out_name,outdf,save_dir),callback=cbk)
            #get_a_ex_a_out(ex_name,ex_path,out_name,out_path,save_dir)
    pool.close()
    pool.join()