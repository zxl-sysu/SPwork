import os
import multiprocessing as mtp
from tqdm import tqdm as tqdm
import shutil
import time
import pandas as pd



def run_cmd(cmd):
    print(cmd)
    os.system(cmd)
    
if __name__=='__main__':

    exposure_path='./match_ex_out'
    exposure_list=os.listdir(exposure_path)
    
    
    for base_name in exposure_list:
        base_path='%s/%s' %(exposure_path,base_name)
        target_list=os.listdir(base_path)
        pool=mtp.Pool()
        for target_name in target_list:

            path='%s/%s/ex_out_sig.txt' % (base_path,target_name)
            
            
            save_dir_forward='foward/%s/%s' %(base_name,target_name)
            if os.path.exists(save_dir_forward):
                shutil.rmtree(save_dir_forward)
                print('removing...... %s' %save_dir_forward)
                time.sleep(2)
            os.makedirs(save_dir_forward,exist_ok=True)
        
            
            cmd='Rscript run_a_MR_forward2.R "%s" "%s" "%s" "%s" "%s" '%(save_dir_forward ,path,"NOlift",base_name,target_name)

            pool.apply_async(run_cmd,args=(cmd,))
    
        pool.close()
        pool.join()

    print('+++++++++++++++++++++++++++++++DONE+++++++++++++++++++++++++++++++')
