import os 
import re
import pandas as pd
import sys

PRS_dir= '.'
disease_name='SPd18' ##改成相同的表型名字

os.chdir(PRS_dir)
print('PRS dir = %s' % PRS_dir)
print('changing working dir to %s'% os.getcwd() )

trait='%s' % (disease_name)

clps=[ i for i in os.listdir() if i.startswith('clump_') ]

pts=[ os.path.join(i,j)+ '/' for i in clps for j in os.listdir(i) ]

print(pts)
gwas=pd.read_csv('./qc_base.txt'  ,sep=' ')
print("%s READ__DONE=== " % PRS_dir )
preffix=trait

for pt in pts:
    vars=[ os.path.join(pt,i) for i in os.listdir(pt) if i.endswith('vars')]
    #print(vars)

    if vars:
        vars.sort(key=lambda x: int(re.findall(r'\d+',x)[0]))
        snp_list=[]
        for i in vars:
            with open(i,'r') as f:
                snp=f.readlines()
                snp_list.append(snp)
        snp_list2=[j.strip('\n') for i in snp_list for j in i] 

        used_snps=gwas[gwas.rsID.isin(snp_list2)]

        used_snps.to_excel('%s/used_gwas.xlsx' %pt  ,index=False)

        print(pt, 'done   \n')
if os.path.exists('all_PRS_result'):
    os.system('rm -r ./all_PRS_result')
if os.path.exists('all_PRS_%s.tar.gz' % preffix):
    os.system('rm-r all_PRS_%s.tar.gz'% preffix)
os.system('mkdir -p all_PRS_result')
os.system("/bin/bash -c 'cp --parents -rf clump_*/*/{sum_score*,used*,hst50*,KS*} all*t/' ") #直接这样写会报错 问题在{} 要加上/bin/bash -c
os.system('tar -zcvf all_PRS_%s.tar.gz all*t' % preffix)
#mkdir -p all_PRS_result
#cp --parents -rf clump_*/*/{sum_score*,used*} all*t/
#tar -zcvf all_PRS.tar.gz all*t

## python multi_clump_PRS.py && python multicp.py