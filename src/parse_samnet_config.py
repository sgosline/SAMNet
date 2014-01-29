'''
Additional file that takes SAMNet arguments from a configuration file instead of the 
command line, then creates a run.sh script to run the appropriate python commands
'''

__author__='Sara JC Gosline'
__email__='sgosline@mit.edu'

import re,sys,os
from ConfigParser import ConfigParser
from optparse import OptionParser


def main():
	
	parser=OptionParser()
	parser.add_option('--conf',dest='conf',type='string',help="Name of configuration file to parse")
#	parser.add_option('--outputdir',dest='outputdir')

	opts,args=parser.parse_args()

        config=ConfigParser()
	config.read(opts.conf)
        print 'Reading file '+opts.conf
        ##these are the samnet parameters
	gamma=config.get('Samnet','gamma')
	tfa=config.get('Samnet','tfa')
	ppi=config.get('Samnet','ppi')
	protweights=config.get('Samnet','protweights')
	mrna=config.get('Samnet','mrna')

	##these are the parameters necessary to run the script properly
	outputdir=config.get('Server','outputdir')
	pythondir=config.get('Server','pythondir')

	davidthresh=config.get('David','threshold')
	

	pythonstr='python '+pythondir+'/samnet.py --PPI='+ppi+' --gamma='+gamma+' --proteinWeightsByComm='+protweights+' --traByComm='+mrna+' --doMCF --tfmrna='+tfa+' --output='+outputdir+'/results/samnetout'
	##now write scrip
        davidstr='python '+pythondir+'/go_enrichment/DAVID.py --outputDir='+outputdir+'/results --listName=samnet --thresholdValue='+davidthresh+' '+outputdir+'/results/samnetoutmultiComm_edge_commodity.eda'

	#os.system('mkdir '+optsoutputdir+'/sh')
	fop=open(outputdir+'/run.sh','w')
	fop.write('source /opt/python2.7.5/env.sh\n')
	fop.write('source /opt/ampl/env.sh\n')
	fop.write(pythonstr+'\n')
	fop.write(davidstr+'\n')

	fop.close()

if __name__=='__main__':
	main()
