'''
This is a script that runs SAMNet on human EMT data

Copyright (c) 2012 Sara JC Gosline, Sarah J Spencer
sgosline@mit.edu, sjspence@mit.edu


'''
__author__='Sara JC Gosline'
__author__='Sarah J Spencer'

__email__='sgosline@mit.edu'
__email__='sjspence@mit.edu'


from optparse import OptionParser
import os,re
import itertools


# EMT data from Thomson et al. 2010.
#initially added by sjs
#updated by sjcg to include running in reverse
def runEMTdata(gamma,cap,doMCF=False,upanddown=False,uselowthresh=False,bp_upstream='5000'):
    '''
    Runs EMT data. if number_to_run >1, then it will run pairs of data sets or 
    all three if doMCF is set then will try mcf version
    '''
    tfa_files={'1000':'9606mitab.wgEncodeOpenChromDnaseA549_1000_upstream.tfa','3000':'9606mitab.wgEncodeOpenChromDnaseA549_3000_upstream.tfa','5000':'9606mitab.wgEncodeOpenChromDnaseA549_5000_upstream.tfa'}
    
    chem_list=['fixed','snail','tgfb','zeb1']
    
    #use expr_files with updated TCF7L2 gene name
    expr_files={'fixed':'exp/HC_exp.txt','snail':'exp/S_exp.txt','tgfb':'exp/T_exp.txt','zeb1':'exp/Z_exp.txt'}

    phen_files={'fixed':'phos/HC_pY_icrogid.phen','snail':'phos/S_pY_icrogid.phen','tgfb':'phos/T_pY_icrogid.phen','zeb1':'phos/Z_pY_icrogid.phen'}
    
    os.system('mkdir emtOutput')

    if uselowthresh: # use lower threshold (takes MUCH longer)
        ppinet='9606.mitab.01192011.uniq_miscore-localirefindex3-20110831.digraph_thresh0.3.pkl'
    else:
        ppinet='9606.mitab.01192011.uniq_miscore-localirefindex3-20110831.digraph_thresh0.5.pkl'
    os_strs=[]
    
    ch=','.join(chem_list)
    exprs=','.join([expr_files[f] for f in chem_list])
    phens=','.join([phen_files[f] for f in chem_list])
    
    output_str='emtOutput/'+re.sub(',','_',ch)+'_network_gamma'+gamma+'_cap_'+cap.split('.')[1]
    
    if uselowthresh:
        output_str+='_Psiquic_3thresh'
    else:
        output_str+='_Psiquic_5thresh'

    exec_str='python ../../src/samnet.py --solver=cplexamp --proteinWeights='+phens+' --tra='+exprs+' --tfmrna='+tfa_files[bp_upstream]+' --updateIds=humaniref --PPI='+ppinet+' --gamma='+gamma+' --treatmentNames='+ch+' --cap='+cap
        
    if doMCF:
        exec_str+=' --doMCF'
            
    if upanddown:
        ostr=output_str+'_upOnly'
        os.system(exec_str+' --output='+ostr+' --upordown=up')
        os_strs.append(ostr)

        ostr=output_str+'_downOnly'
        os.system(exec_str+' --output='+ostr+' --upordown=down')
        os_strs.append(ostr)
        
    else:
        os.system(exec_str+' --output='+output_str)
        os_strs.append(output_str)
        
    return os_strs


def emt_stats():

    chem_list=['fixed','snail','tgfb','zeb1']
    #use expr_files with updated TCF7L2 gene name
    expr_files={'fixed':'exp/HC_exp.txt','snail':'exp/S_exp.txt','tgfb':'exp/T_exp.txt','zeb1':'exp/Z_exp.txt'}

    phen_files={'fixed':'phos/HC_pY_icrogid.phen','snail':'phos/S_pY_icrogid.phen','tgfb':'phos/T_pY_icrogid.phen','zeb1':'phos/Z_pY_icrogid.phen'}

    file_out=open('emt_overlap_stats.xls','w')
    file_out.write('Commodity\tPhosphoprotein hits\tDiff Ex Genes\tOverlap\n')

    for metal in chem_list:
        row=metal
        exp=[arr.split('\t')[0].strip() for arr in open(expr_files[metal],'r').readlines()]
        phen=[arr.split('\t')[0].strip() for arr in open(re.sub('_icrogid','',phen_files[metal]),'r').readlines()]
        print 'Got '+str(len(exp))+' expression values and '+str(len(phen))+' phen values for '+metal
        overlap=[p for p in phen if p in exp]
        row+='\t'+str(len(phen))+'\t'+str(len(exp))+'\t'+str(len(overlap))+'\n'
        file_out.write(row)
    file_out.close()    


   
def main():
    parser=OptionParser()
    parser.add_option('--gamma',type='string',dest='gamma',help='Gamma value to use. DEFAULT=14',default='14')

    parser.add_option('--dataSet',type='string',dest='dataset',help='OPTIONAL: Set to emt_stats to analyze dataset statistics',default='emt')
    
    parser.add_option('--bpUpstream',type='string',dest='bp_upstream',help='Set to distance upstream. Choice of 1000,3000 or 5000. Can also add comma-delimited list. Default: 5000', default='5000')
    
    parser.add_option('--cap',type='string',dest='cap',help='Set to cap for ppi network. Default is .99',default='.99')

    parser.add_option('--doRN',action='store_true',dest='mcf',help='Set to run original ResponseNet. DEFAULT is False',default=False)
    parser.add_option('--upanddown',action='store_true',dest='upanddown',help='Set to run up and down-regulated genes individually',default=False)
    parser.add_option('--useLowThresh',action='store_true',dest='usestring',help='Set if you want to use include lower confidence PSIQUIC interactions',default=False)



    (opts,args)=parser.parse_args()
    if opts.dataset=='emt':
        for bp in opts.bp_upstream.split(','):
            output_files=runEMTdata(gamma=opts.gamma,cap=opts.cap,doMCF=not opts.mcf,upanddown=opts.upanddown,uselowthresh=opts.usestring,bp_upstream=bp)
    elif opts.dataset=='emt_stats':
        output_files=emt_stats()
        
if __name__=='__main__':
    main()
