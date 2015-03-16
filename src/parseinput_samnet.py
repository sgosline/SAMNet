'''
parseinput_samnet.py
SAMNet module handles pre-processing of SAMNet files

Copyright (c) 2012 Sara JC Gosline
sgosline@mit.edu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

'''

import re,math,networkx,collections,pickle

def get_weights_phen_source(phendatadict):   

    weights_source = collections.defaultdict(dict)
    for treatment in phendatadict.keys(): #changed for multi-source
        phendatalist=phendatadict[treatment]
        for item in phendatalist:
        #weight associated with the protein(from phenotypic data) and the source
        #this weight is given in somephenfile (in the column following the protein)
            fields = item.strip('\r\n').split('\t')
            protein = fields[0].strip()
            weight =math.fabs(float(fields[1].strip()))
            if weight==0.0:
                weight=0.000001
            if weight==1.0:
                weight=0.999999
            weights_source[treatment][protein] = weight
    
    #keys of this dictionary will be all proteins from the genes of the phenotypic dataset
    return weights_source

def by_comm_into_one_dict(argument_file_containing_comm,doUpper=False):
    '''
    Takes multiple file arguments but parses them into separate elements of a dictionary
    New function for multi-source
    '''
    argument = argument_file_containing_comm
        
    lines=[line.strip() for line in open(argument)]
    alldata=dict()
    treats=set()
    mrna=set()
    for l in lines:
        arr=l.strip().split()
        m=arr[1].strip()
        if doUpper:
            m=m.upper()
        treats.add(arr[0])
        if arr[0]+'_treatment' not in alldata:
            alldata[arr[0]+'_treatment']=[]
        alldata[arr[0]+'_treatment'].append(m+'\t'+arr[2])
        mrna.add(m)

    #print ','.join(mrna)
    return alldata,[a for a in treats],[a for a in mrna]

def get_direct_target_weights(targ_dict):
    '''
    Takes as input a dictionary of tab-delimited files (opened as list) with mRNA targets and weights
    for a single miRNA
    '''

    targ_dict=collections.defaultdict(dict)#changed for multi source
    for treatment in targ_dict.keys():
        targs=targ_dict[treatment]
        for line in targs:
            [mrna,weight]=re.split('\t',line.strip())
            w=float(weight)
            if(w==0.0):
                w=0.000001
            if w==1.0:
                w=0.999999
            targ_dict[treatment][mrna+'mrna']=math.fabs(w)

    return targ_dict

def multiple_args_into_one_list(arguments_separated_by_comma):
        arguments = arguments_separated_by_comma
        #make a list of all the arguments given
        list_of_datasets = arguments.split(",")
        #initialize the list that will hold data from all datasets

        alldata = []
        if len(list_of_datasets)==0 or list_of_datasets=='':
            return alldata
        #go through the list of datasets
        for dataset in list_of_datasets:
            if dataset=='':
                continue
            onedataset = open (dataset, 'r')
            #put each line of the dataset into a list that will contain data from all datasets

            lines = onedataset.readlines()
            for line in lines:
                if '@' not in line:
                    alldata.append(line)

        return alldata
    
def multiple_args_into_one_dict(arguments_separated_by_comma,treatmentNames=''):
    '''
    Takes multiple file arguments but parses them into separate elements of a dictionary
    New function for multi-source
    '''
    arguments = arguments_separated_by_comma
        #make a list of all the arguments given
    list_of_datasets = arguments.split(",")
    alldata = collections.defaultdict(list)

    if len(list_of_datasets)==0 or list_of_datasets=='':
        return alldata
        #initialize the list that will hold data from all datasets
        
    if(len(list_of_datasets)>1):
        treatlist=['treatment'+str(count) for count in range(len(list_of_datasets))]
    else:
        treatlist=['treatment']

    allt=treatmentNames.split(',')
    if(len(allt)!=len(list_of_datasets)):
        print "List of treatments does not equal list of files, using default naming"
            
    elif len(allt[0])>0:
        for i in range(len(allt)):
            treatlist[i]=allt[i]+'_treatment'



        #go through the list of datasets
    for dataset in range(len(list_of_datasets)):
        if list_of_datasets[dataset]=='':
            continue
        onedataset = open (list_of_datasets[dataset], 'r')
            #put each line of the dataset into a list that will contain data from all datasets

        lines = onedataset.readlines()
        for line in lines:
            if '@' not in line:
                alldata[treatlist[dataset]].append(line)

    return alldata

def filter_ppi(ppi_network,list_of_prots,updateIds=''):
    '''
    Filters protein interaction network by list of proteins for tissue-specific analysis
    '''
    prots=[i.strip() for i in list_of_prots.readlines()]
#    print prots[1:10]
    print "Filtering network of with",str(ppi_network.number_of_nodes()),"nodes and",str(ppi_network.number_of_edges()),"edges with list of",str(len(prots)),"proteins"

    for p in ppi_network.nodes():
        if(p not in prots):
            #ppi_network.delete_node(p)
            ppi_network.remove_node(p)
    print "Network now has",str(ppi_network.number_of_nodes()),'nodes and',str(ppi_network.number_of_edges()),'edges'
    return ppi_network

'''
These two functions normalize responseNet dictionaries or networks
'''

def renormalizeNetworkweights(ppi,type='ecdf'):
    '''
    Takes networkx encoded graph and replaces edge weights with those that are normalized to fall between
    0 and 1
    '''
    all_weights=[ppi.get_edge_data(x[0],x[1])['weight'] for x in ppi.edges()]
    #for each unique weight create a dictionary lookup to replace weight with normalized weight
    wdict={}
    for w in all_weights:
        if w not in wdict.keys():
            wdict[w]=float(len([y for y in all_weights if y<=w]))/float(len(all_weights))
    #iterate through and replace edgeweights
    for x in ppi.edges():
        ppi.get_edge_data(x[0],x[1])['weight']=wdict[ppi.get_edge_data(x[0],x[1])['weight']]
    return ppi



def renormalizeDictionaryweights(ppi,type='ecdf'):
    '''
    Takes networkx encoded graph and replaces edge weights with those that are normalized to fall between
    0 and 1
    '''
    all_weights=ppi.values()
   # print all_weights[0:5]
    #for each unique weight create a dictionary lookup to replace weight with normalized weight
    wdict={}
    for w in all_weights:
        sw=str(w)
        if sw not in wdict.keys():
            if type=='ecdf':
                wdict[sw]=float(len([y for y in all_weights if y<=w]))/float(len(all_weights))
            elif type=='probability': ##this assumes they sum to 1
                wdict[sw]=float(w)/float(sum(all_weights))
            elif type=='fraction': ##fraction of max
                wdict[sw]=float(w)/float(max(all_weights))
   # print wdict
    #iterate through and replace edgeweights
    for x in ppi.keys():
        if str(ppi[x]) in wdict.keys():
            ppi[x]=wdict[str(ppi[x])]

    return ppi


def get_ppi_network(networkfile,doUpper=False):
    '''
    First figure out if the file is a pickle or a text file, then make into network if text
    networkfile: text or pickle file
    doupper: make uppercase if try and in text file
    '''
    #figure out file format
    ext=networkfile.split('.')[-1]

#        print 'Determined '+networkfile+' is PKL, loading...'
    try:
        ppinet=pickle.load(open(networkfile,'rU'))
        ##now test for object type
        if type(ppinet) is networkx.Graph or type(ppinet) is networkx.DiGraph:
            print 'Unpickled is networkX, returning'
        else:
            print 'Pickled object unknown, returning empty digraph'
            ppinet=networkx.DiGraph()

    except:
        print networkfile+' is Not pickle, assuming text'
        ppinet=networkx.DiGraph()
        for row in open(networkfile,'rU').readlines():
            arr=row.strip().split()
            if len(arr)!=3:
                continue
            prot1=arr[0]
            prot2=arr[1]
            if doUpper:
                prot1=prot1.upper()
                prot2=prot2.upper()
            weight=float(arr[2])
            if weight==1.0:
                weight=0.999999
            ppinet.add_edge(prot1,prot2,{'weight':weight})
            ppinet.add_edge(prot2,prot1,{'weight':weight})
        
    

    print "Returning "+str(ppinet.number_of_edges())+' PPI interaction scores'

    return ppinet
