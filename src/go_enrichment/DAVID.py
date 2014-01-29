
'''
Uses web service to do DAVID analysis for cuff diff and MISO parsed output

'''


import re,sys,os
from optparse import OptionParser
from collections import defaultdict


# import ChartReport

##use uniprot dict
uniprot_map='/usr/local/apachedocs/htdocs/samnetwebhome/sample_data/uniprot-human.tab'
up_dict=defaultdict(list)
for row in open(uniprot_map,'r').readlines()[1:]:
	arr=row.strip().split('\t')
	allgenes=arr[3].split()
	for a in allgenes:
		up_dict[a]=arr[0]
		
		
#from tests import *
#from suds import *
#import logging
def getGenesFromCuffDiffFile(filename,pvalue=0.01):
	'''
	gets gene list from cuffdiff file
	'''
	genes=[]

	genes=[arr.split('\t')[0] for arr in open(filename,'r').readlines()[1:] if float(arr.split('\t')[11])<pvalue]

	mapped=[]
	missed=[]
	for g in genes:
		if g in up_dict.keys():
			mapped.append(up_dict[g])
		else:
			missed.append(g)
			
	print 'Got '+str(len(mapped))+' up ids from '+str(len(genes))+' and missed '+','.join(missed)
	return mapped


def getGenesFromMisoParsedFile(filename):
	'''
	gets gene list from miso parsed file
	'''
	genes=[]
	allrows=open(filename,'r').readlines()
	cols=allrows[0].strip().split('\t')
	geneind=cols.index('Gene_ID')

	genes=[re.sub(',','\n',row.strip().split('\t')[geneind]) for row in allrows[1:]]
	
	return genes


def getGenesFromAPclusterOutput(filename):
	'''
	gets list of gene lists from ap cluster file
	'''
	genes=defaultdict(list)
	allrows=open(filename,'r').readlines()
	cols=allrows[0].strip().split('\t')
	geneind=cols.index('gene name')
	clustind=cols.index('cluster number')
	
	for row in allrows[1:]:
		arr=row.strip().split('\t')
		genes[arr[clustind]].append(arr[geneind])

	mapped=defaultdict(list)
	bg=set()
	missed=[]
	for clust in genes.keys():
		for g in genes[clust]:
			if g in up_dict.keys():
				mapped[clust].append(up_dict[g])
				bg.add(up_dict[g])
			else:
				missed.append(g)
			
	print 'Got '+str(len(mapped.values()))+' up ids from '+str(len(genes.values()))+' and missed '+','.join(missed)+' with a bg size of '+str(len(bg))
	return mapped,bg

def getGenesFromCytoscapeSIFFile(filename):
	'''
	gets gene list from Cytoscape sif file
	'''

	genes=defaultdict(list)
	allrows=open(filename,'r').readlines()
	nonsourcesink=[]
	sourcesink=[]
	for row in allrows:
	    arr=row.strip().split('\t')
	    if ('S1' not in arr) and ('T1' not in arr):
		nonsourcesink.append(arr)
		genes[arr[1]].append(arr[0])
		genes[arr[1]].append(arr[2])
	    else:
		sourcesink.append(arr)

	mapped=defaultdict(list)
	bg=set()
	missed=[]
	for mir in genes.keys():
		nonUniqueGenes=genes[mir]
		nonUniqueGenesSet=set(nonUniqueGenes)
		uniqueGenes=list(nonUniqueGenesSet)
		genes[mir]=uniqueGenes
	
		for g in genes[mir]:
		        if g in up_dict.keys():
		                mapped[mir].append(up_dict[g])
		                bg.add(up_dict[g])
		        else:
		                missed.append(g)
		                
	missed = set(missed)
	missed = list(missed)


	print 'Got '+str(len(mapped.values()))+' up ids from '+str(len(genes.values()))+' and missed '+','.join(missed)+' with a bg size of '+str(len(bg))
	return mapped,bg        
	
		
def getGenesFromCytoscapeEDAFile(filename):
	'''
	gets gene list from Cytoscape eda file
	'''

	genes=defaultdict(list)
	allrows=open(filename,'r').readlines()[1:]
	nonsourcesink=[]
	sourcesink=[]
	for row in allrows:
	    arr=row.strip().split(' ')
	    if ('S1' not in arr) and ('T1' not in arr):
		nonsourcesink.append(arr)
                comm=arr[1].strip('(').strip(')')
		genes[comm].append(arr[0])
		genes[comm].append(arr[2])
	    else:
		sourcesink.append(arr)

	mapped=defaultdict(list)
	bg=set()
	missed=[]
	for mir in genes.keys():
		nonUniqueGenes=genes[mir]
		nonUniqueGenesSet=set(nonUniqueGenes)
		uniqueGenes=list(nonUniqueGenesSet)
		genes[mir]=uniqueGenes
	
		for g in genes[mir]:
		        if g in up_dict.keys():
		                mapped[mir].append(up_dict[g])
		                bg.add(up_dict[g])
		        else:
		                missed.append(g)
		                
	missed = set(missed)
	missed = list(missed)

	print 'Got '+str(len(mapped.values()))+' up ids from '+str(len(genes.values()))+' and missed '+','.join(missed)+' with a bg size of '+str(len(bg))
	return mapped,bg        	

def callDavidCommand(genelist,listname,daviddir):
	'''
	do i need this?
	'''
	cmd='python '+daviddir+'/DAVIDWebService_Client.py --listName='+listname+' '+','.join(genelist)
	print cmd
	os.system(cmd)


def doEnrichment(genefile,listname,bgfile,bgname='',idType='OFFICIAL_GENE_SYMBOL',daviddir=''):
	categories='GOTERM_BP_FAT,GOTERM_MF_FAT,GOTERM_CC_FAT,BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'
	if bgfile=='':
		cmd='python '+daviddir+'/ChartReport.py --idType='+idType+' --listName='+listname+' '+genefile
	else:
		cmd='python '+daviddir+'/ChartReport.py --idType='+idType+' --listName='+listname+' '+genefile+' '+bgfile
	print cmd
	os.system(cmd)


 #   ChartReport.DAVIDenrich(genefile,idType=idType,bgF=bgfile,bgName=bgname,listName=listname,category=categories)
	
def main():
	parser=OptionParser()
	parser.add_option('--listName',dest='listName',type='string',default='',help='Name of list')
	parser.add_option('--listType',dest='listType',type='string',default='CYTOEDA',help='MISO or cuffdiff file to parse?')
	parser.add_option('--threshold',dest='threshold',type='string',default='Benjamini',help='Threshold for significant DAVID terms: Count or % or Pvalue or List Total or Pop Hits or Pop Total or Fold Enrichment or Bonferroni or Benjamini or FDR')
	parser.add_option('--thresholdValue',dest='thresholdValue',type='float',default='0.1',help='Threshold less than value'
	)
	parser.add_option('--outputDir',dest='outputDir',type='string',default='./',help='Directory for output chartReport and xls file')
					  
	opts,args=parser.parse_args()
        daviddir=os.path.dirname(sys.argv[0])
	print daviddir
	
	## Check if output directory exists ##
	try:
    		os.makedirs(opts.outputDir)
	except OSError:
    		pass
        
	if opts.listType=='MISO':
		genes=getGenesFromMisoParsedFile(args[0])
		idtype='ENSEMBL_GENE_ID'
		fname=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'_genes.txt'
		open(fname,'w').writelines([g+'\n' for g in genes])
		flist=[fname]
		
	elif opts.listType=='CUFFDIFF':
		genes=getGenesFromCuffDiffFile(args[0])
		idtype='UNIPROT_ACCESSION'
		fname=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'_genes.txt'
		open(fname,'w').writelines([g+'\n' for g in genes])
		flist=[fname]
		
	elif opts.listType=='APCLUSTER':
		genes,bg=getGenesFromAPclusterOutput(args[0])
		idtype='UNIPROT_ACCESSION'
		flist=[]
		for clustid in genes.keys():
			fname=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'clust_'+clustid+'_genes.txt'
			flist.append(fname)
			open(fname,'w').writelines([g+'\n' for g in genes[clustid]])

	elif opts.listType=='CYTOSIF':
		genes, bg =getGenesFromCytoscapeSIFFile(args[0])
		idtype='UNIPROT_ACCESSION'
		flist=[]
		for mir in genes.keys():
			fname=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'_'+mir+'_genes.txt'
			flist.append(fname)
			open(fname,'w').writelines([g+'\n' for g in genes[mir]])

	elif opts.listType=='CYTOEDA':
		genes, bg =getGenesFromCytoscapeEDAFile(args[0])
		idtype='UNIPROT_ACCESSION'
		flist=[]
		for mir in genes.keys():
			fname=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'_'+mir+'_genes.txt'
			flist.append(fname)
			open(fname,'w').writelines([g+'\n' for g in genes[mir]])
	
	
	if len(args)>1:
		if opts.listType=='MISO':
			bggenes=getGenesFromMisoParsedFile(args[1])

		else:
			bggenes=getGenesFromCuffDiffFile(args[1])

		bgfile=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'_BG_genes.txt'
		bgname='background'

	elif opts.listType=='APCLUSTER': #we have defacto background
		bgfile=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'_BG_genes.txt'
		bgname='all_clustered_Genes'
		open(bgfile,'w').writelines([b+'\n' for b in bg])
	
	elif opts.listType=='CYTOSIF':	#we have defacto background
		bgfile=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'_BG_genes.txt'
		bgname='all_clustered_Genes'
		open(bgfile,'w').writelines([b+'\n' for b in bg])

	elif opts.listType=='CYTOEDA':	#we have defacto background
		bgfile=opts.outputDir+'/'+opts.listName+'_'+opts.listType+'_BG_genes.txt'
		bgname='all_clustered_Genes'
		open(bgfile,'w').writelines([b+'\n' for b in bg])


	else:
		bggenes=[]
		bgfile=''
		bgname=''
  
	if(bgfile==''):
		for fname in flist:
			doEnrichment(os.path.realpath(fname),opts.listName,bgfile,bgname,idType=idtype,daviddir=daviddir)
	else:
		for fname in flist:
			doEnrichment(os.path.realpath(fname),opts.listName,os.path.realpath(bgfile),bgname,idType=idtype,daviddir=daviddir)
			
	## Collect Sig DAVID Terms ##
	davidfiles=[f for f in os.listdir(opts.outputDir) if 'chartReport' in f]
	print davidfiles,"\n\n"
	print 'Got '+str(len(davidfiles))+' files with DAVID Results'

	#correct='Benjamini'
	correct=opts.threshold
	
	final_file=open(opts.listName+'_'+opts.listType+'_'+opts.threshold+'_'+str(opts.thresholdValue)+'_sigDavidTerms.xls','w')
	final_file.write('origFile\tCategory\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\t Bonferroni\t Benjamini\tFDR\n')
	
	for df in davidfiles:
		allrows=open(opts.outputDir+'/'+df,'r').readlines()
		colnames=allrows[0].strip().split('\t')
		if correct not in colnames:
			print colnames
			continue
		bh=colnames.index(correct)##can switch between other indices
		for row in allrows[1:]:
			if float(row.strip().split('\t')[bh])<opts.thresholdValue:
				final_file.write(df.split('.txt.')[0]+'\t'+row)
	final_file.close()
	
	
if __name__=='__main__':
	main()
