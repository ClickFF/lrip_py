import shutil
import re
import os


def prep_gb_decomp(ligs, gb_min_path, dec_in_path, comp_info, min_end_idx):
    
    current_path = os.getcwd()
    #gb_min_path = current_path + "/GB_MIN_%s"%num_ite
    gb_dec_path = current_path
    #dec_in_path = dec_in_path + "/gb_decomp"
    recsrt = comp_info['recsrt']
    recend = comp_info['recend']
    ligsrt = comp_info['ligsrt']
    ligend = comp_info['ligend']
    
    for lig in ligs['id']:
        if os.path.isdir('%s/%s'%(gb_dec_path,lig)):
            pass
        else:
            os.mkdir('%s/%s'%(gb_dec_path,lig))
        if os.path.isdir('%s/%s/minout'%(gb_dec_path,lig)):
            pass
        else:
            os.mkdir('%s/%s/minout'%(gb_dec_path,lig))
        shutil.copy('%s/%s/prmtop'%(gb_min_path,lig), '%s/%s/prmtop'%(gb_dec_path,lig))
        shutil.copy('%s/%s/min%s.rst'%(gb_min_path,lig, min_end_idx), '%s/%s/inpcrd'%(gb_dec_path,lig)) # This mechanism does not support min stages greater than 5 since the RUN.template file is hardcoded to min5.rst.
        with open(dec_in_path, 'r') as min_in:
            min_in_contents = min_in.read()
        
        min_in_contents = re.sub('%RRES%', 'RRES %s  %s'%(recsrt,recend), min_in_contents)
        min_in_contents = re.sub('%LRES%', 'LRES %s  %s'%(ligsrt,ligend), min_in_contents)
        min_in_contents = re.sub('%RES%', 'RES %s  %s'%(recsrt,ligend), min_in_contents)
        with open('%s/%s/min.in'%(gb_dec_path,lig), 'w') as min_in:
            min_in.write(min_in_contents)

if __name__ == '__main__':
    ligs={'d4r_3': 0}
    if_ite = False
    num_ite = 0
    dec_in_path = '/data6/tan77/gpcr/docking/test/docs'
    comp_info = {'recsrt':1,
                 'recend':373,
                 'ligsrt':374,
                 'ligend':374}
    prep_gb_decomp(ligs, if_ite, num_ite, dec_in_path, comp_info)

        


    

