'''
This is a script that runs SAMNet code using the Yeast metal dataset from Jin et al.
Copyright (c) 2012 Sara JC Gosline
sgosline@mit.edu
'''
__author__='Sara JC Gosline'
__email__='sgosline@mit.edu'



from optparse import OptionParser
import os,re
       

def runMetalData(gamma,cap,doMCF=False,do_unif=False,do_ecdf=False,upanddown=False):
    '''
    Runs metal data.  If number_to_run is greater than 1, will put that many metals side by side
    If doMCF is true will run MCF version
    '''

    metal_list=['arsenic','cadmium','chromium','copper','mercury','silver','zinc']
    phen_conc={'arsenic':1250,'cadmium':25,'chromium':900,'copper':7000,'mercury':49,'silver':20,'zinc':2000}
    expr_conc={'arsenic':1250,'cadmium':25,'chromium':1700,'copper':9000,'mercury':47,'silver':20,'zinc':2000}
    #expr_conc={'arsenic':1250,'cadmium':25,'chromium':400,'copper':9000,'mercury':47,'silver':20,'zinc':2000}
    os.system('mkdir metalOutput')    


    os_strs=[]

    #print c
    mets=','.join(metal_list)
    exps=','.join(['./exp/'+metal+'_'+str(expr_conc[metal])+'jinEtAlExpscores'+'_001'+'.txt' for metal in metal_list])
        
    #output string prefix
    output_str='metalOutput/'+re.sub(',','_',mets)+'Yeast_'+gamma+'p01expr'+'_cap'+cap.split('.')[1]
    
    output_str+='_Stringv9_6thresh'

    #            output_str+='_origPPI_withkinase'
    
    #phenotypic files. each set is normalized slightly differently
    phens=','.join(['./phen/'+metal+'_'+str(phen_conc[metal])+'jinEtAl_EC50_PhenScores_5thresh.phen' for metal in metal_list])
    if do_unif:
        phens=','.join(['./phen/'+metal+'_'+str(phen_conc[metal])+'jinEtAl_EC50_PhenScores_unif_scores_5_thresh.phen' for metal in metal_list])
        output_str+='_unif_phens'
    elif do_ecdf:
        phens=','.join(['./phen/'+metal+'_'+str(phen_conc[metal])+'jinEtAl_EC50_PhenScores_5thresh_ecdf_norm.phen' for metal in metal_list])
        output_str+='_ecdf_phens'

    ##we include two protein-protein interaction networks
    ppinet='4932.string.9.0.ExpScores6thresh.pkl'
            
    if not upanddown:
        exec_str='python ../../src/samnet.py --proteinWeights='+phens+' --tra='+exps+' --PPI='+ppinet+' --tfmrna=allWeighted_tf2gene.tfa --gamma='+gamma+' --output='+output_str+' --updateIds=Yeast --cap='+cap+' --treatmentNames='+mets+' --debug --solver=cplexamp'
        if doMCF:
            exec_str+=' --doMCF'
        os.system(exec_str)

    else:##run up only and down only 
        exec_str='python ../../src/samnet.py --proteinWeights='+phens+' --tra='+exps+' --PPI='+ppinet+' --tfmrna=allWeighted_tf2gene.tfa --gamma='+gamma+' --output='+output_str+'_upOnly'+' --updateIds=Yeast --cap='+cap+' --treatmentNames='+mets+' --upordown=up --debug --solver=cplexamp'
        os_strs.append(output_str+'_upOnly')
        
        if doMCF:
            exec_str+=' --doMCF'
        os.system(exec_str)
        
        exec_str='python ../../src/samnet.py --proteinWeights='+phens+' --tra='+exps+' --PPI='+ppinet+' --tfmrna=.allWeighted_tf2gene.tfa --gamma='+gamma+' --output='+output_str+'_downOnly'+' --updateIds=Yeast --cap='+cap+' --treatmentNames='+mets+' --upordown=down --debug --solver=cplexamp'
            
        if doMCF:
            exec_str+=' --doMCF'
        os.system(exec_str)
        os_strs.append(output_str+'_downOnly')

    return os_strs



def metal_stats():
    metal_list=['arsenic','cadmium','chromium','copper','mercury','silver','zinc']
    phen_conc={'arsenic':1250,'cadmium':25,'chromium':900,'copper':7000,'mercury':49,'silver':20,'zinc':2000}
    expr_conc={'arsenic':1250,'cadmium':25,'chromium':1700,'copper':9000,'mercury':47,'silver':20,'zinc':2000}
    
#    mets=','.metal_list
    file_out=open('metal_overlap_stats.xls','w')
    file_out.write('Commodity\tPhenotypic hits\tDiff Ex Genes\tOverlap\n')

    for metal in metal_list:
        row=metal
        exp=[arr.split('\t')[0].strip() for arr in open('./exp/'+metal+'_'+str(expr_conc[metal])+'jinEtAlExpscores'+'_001'+'.txt','r').readlines()]
        phen=[arr.split('\t')[0].strip() for arr in open('./phen/'+metal+'_'+str(phen_conc[metal])+'jinEtAl_EC50_PhenScores_5thresh.phen','r').readlines()]
        
        print 'Got '+str(len(exp))+' expression values and '+str(len(phen))+' phen values for '+metal
        overlap=[p for p in phen if p in exp]
        row+='\t'+str(len(phen))+'\t'+str(len(exp))+'\t'+str(len(overlap))+'\n'
        file_out.write(row)
        
    file_out.close()    


   
def main():
    parser=OptionParser()
    parser.add_option('--gamma',type='string',dest='gamma',help='Gamma value to use. DEFAULT=15',default='15')
    parser.add_option('--dataSet',type='string',dest='dataset',help='OPTIONAL: set to metal_stats to just print out dataset statistics.',default='metal')
    parser.add_option('--cap',type='string',dest='cap',help='OPTIONAL: Set to cap for ppi network. Default is .99',default='.99')
    parser.add_option('--doRN',action='store_true',dest='mcf',help='Set to run original ResponseNet. DEFAULT is False',default=False)
    parser.add_option('--upanddown',action='store_true',dest='upanddown',help='Set to run up and down-regulated genes individually. DEFAULT: False',default=False)
    

    (opts,args)=parser.parse_args()

    if opts.dataset=='metal_random':
        output_files=runMetalRandomization()

    if opts.dataset=='metal':
        output_files=runMetalData(opts.gamma,opts.cap,not opts.mcf,upanddown=opts.upanddown)
    elif opts.dataset=='metal_stats':
        output_files=metal_stats()

        


if __name__=='__main__':
    main()

