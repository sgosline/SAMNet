'''
samnet.py
Primary samnet executable

Copyright (c) 2012-2014 Sara JC Gosline
sgosline@mit.edu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

'''
import sys,pickle
from numpy import *
import os,re
##new addition 2013-01-17, replacing os.system with spawn command
import subprocess
import networkx
from optparse import OptionParser
import time
from collections import defaultdict


#This file parses the command line input, executes the ampl command and combines multiple commodities if MCF is not used

#these modules do the rest:
import parseinput_samnet as parseIn##handles input and file parsing and filtering
import writefiles_samnet as wf #handles writing files to be run by ampl
import post_samnet as post ##handles post-processing of ampl results
import tfnetwork_samnet as tfNetwork #handles transcriptional data

fpath=os.path.dirname(os.path.abspath( __file__ ))
print fpath

id_directory=re.sub('src','lib',fpath)
print id_directory

def main():
    #PARSING PARAMETERS
    parser = OptionParser()

    #############################################################
    ###INPUTs to flow network: phenotypic hits/indirect targets
    parser.add_option("--proteinWeightsByComm", type="string", dest="prot_weights_comm", help="Phenotypic response file (if using genetic data). '.txt' format . Must contain column of treatments (commodities), column of *protein* names and column of weights, delimited by tab. For multiple treatments, only one file should be provided.",default='')

    parser.add_option("--proteinWeights", type="string", dest="prot_weights", help="Phenotypic response file (if using genetic data). '.txt' format . Must contain column of *protein* names and column of weights, delimited by tab. If multiple files are provided, they will be treated as different treatments and named following the list in --treatmentNames",default='')
    
#    parser.add_option("--phen", type="string", dest="phenfile", help="DEPRECATED, please use --proteinWeights.  Phenotypic response file (if using genetic data). '.txt' format . Must contain column of *protein* names and column of weights, delimited by tab. If multiple files are provided, they will be treated as different treatments and named following the list in --treatmentNames",default='')
    
    parser.add_option("--mrnaWeights",type="string",dest='direct_weights',help="Weight file to model direct effects on mRNA (e.g. miRNA).  '.txt' format. Must contain a column of mRNA targets (gene names must match --tra argument)  and a column of scores depicting how strong a weight should be placed on that mRNA. Absolute value is taken. If multiple files are provided, they will be treated as different treatments and named following the list in --treatmentNames",default='')

#    parser.add_option('--treatmentNames',type='string',dest='treat_names',help='List of names for each treatment provided in the proteinWeights (and mrnaWeights) flag.  If no list is provided treatments are numbered',default='')

    parser.add_option('--treatmentWeights',type='string',dest='treat_weights',help='List of treatment weights, delimited by comma, corresponding to the protein weight files and treatment names (and mrna weights if included). Default values are 1.0s',default='')

    ##############################################################

    ###THE MIDDLE NETWORK: Protein interaction network, string ids
    parser.add_option("--PPI", type="string",  dest="PPIfile", help="Directed PPI graph (of class XDiGraph from networkx) in pkl form. Edges must contain weight information (node1, node2, {'weight':weightvalue}). Undirected edges should be represented as 2 directed edges between the respective 2 nodes.")# If multiple files will be used, enter these at the commandline, one after the other, separated by commas")

    #############################################################

    ####Protein-DNA network
    parser.add_option("--tfmrna", type="string", dest="tf_mrna_weights", help="Weight file for transcriptional network. '.txt' format. Must contain a column of TF (*protein* name), a column of mRNA regulated by that TF (*gene name*) and a column of weight, separated by tab. If multiple files will be used, enter these at the commandline, one after the other, separated by commas",default='')

    ####Transcriptional response file, at this point only one
    parser.add_option("--tra", type="string",  dest="trafile", help="Transcriptional response file. \'.txt\' format. Must contain a column of differentially expressed genes by *gene name*, column of fold difference, and a column of p-values, delimited by tab. If multiple files will be used, enter these at the commandline, one after the other, separated by commas",default='')

    parser.add_option("--traByComm",type='string',dest='trafile_comm',help="Transcriptional response file. \'.txt\' format. Must contain a column of treatments (commodities), a column of differentially expressed genes by *gene name*, and a column of scores depicting how strong a weight should be placed on that mRNA (or instead, have a column of fold difference and a column of p-values), delimited by tab. For multiple treatments, only one file should be provided.",default='')

    
    ##Use this option to exclude mRNA from network...
#    parser.add_option('--nomrna',type='string',dest='no_mrna',help="Set to 'true' to remove mrna from network and calculate weights of transcription factors independently.  Will print out mrna network but not use it to calculate flow", default='False')
    ###AVOID ALL TF-mRNA interactions, or use some other format
#    parser.add_option("--rawTfData",type='string',default='',dest='raw_tf_data',help="OPTIONAL: set to file location of raw transcription factor weights. This option over-rides --tfmrna and --tra to provide new weights to a set of sink proteins. DEPRACATED: use noTfs flag")

    ################################################################
    

    ##GENERAL PARAMETERS for SAMNet
    parser.add_option("--output", type="string", dest="outputfile", help="Output file name. Directory + output file name without any extension")
    parser.add_option('--doMCF',action='store_true',dest='mcf',help='Set this option to run multiple sources/sinks as multiple commodities.DEFAULT is False',default=False)

    parser.add_option('--makeCombined',action='store_true',dest='make_combined',default=False,help='Set this option to create a merged single commodity flow network as well as the mcf.  The files will have a MERGED_FROM_SINGLE prefix ahead of the \'output\' string and no \'multicomm\' on the end.')
    
#    parser.add_option("--foldchangeorp", type="string", default="foldchange", dest="foldchange", help="OPTIONAL: whether to use fold difference or p-value for weights between mRNA and sink. Write 'foldchange' for fold difference, or 'pvalue' for p-value. Default is 'foldchange'")

#    parser.add_option("--foldtraorp", type="string", default="", dest="foldtra", help="DEPRECATED: please use --foldchangeorp.  OPTIONAL: whether to use fold difference or p-value for weights between mRNA and sink. Write 'foldtra' for fold difference, or 'pvalue' for p-value. Default is ''")

    parser.add_option('--upordown',type='string',default='',dest='upordown',help="OPTIONAL: whether or not to use up- or down-regulated genes only. Defaults to using the absolute value, but use 'up' to use upregulated genes or 'down' to use down-regulated genes")

    parser.add_option("--gamma", type="string", default='8', dest="gamma", help="OPTIONAL:Integer number for gamma. Default is 8")
    parser.add_option("--cap", type="float", default=0.99, dest="cap", help="OPTIONAL:Value for the capping of weights. Default is 0.99")

    parser.add_option("--solver", type="string", dest="solver", help="OPTIONAL: Solver used to generate the flow solution. Write exact name of the solver as it should appear in an ampl command. Default is cplexamp", default="cplexamp")
    
    parser.add_option("--updateIds",type='string',dest='updateIds',help="OPTIONAL: Set to \'mouse\' if you want to map protein ids to mouse gene names, \'human\' if you want to map to human gene names. Underlying assumption is that protein-protein interaction network, TF-DNA interaction and source weights are all in UNIPROT identifies.",default='')

    parser.add_option('--hier-cap',action='store_true',dest='hier_cap',help='OPTIONAL: Set this flag to make capacities decrease by factor of 1/10^x where x is unweighted shortest path to Source.  WARNING: not tested for MCF yet.',default=False)
    parser.add_option('--noTfs',action='store_true',dest='noTfs',default=False,help="OPTIONAL: Set this flag to take sink weights and directly connect them to protein interaction network without using a TF-DNA interaction network as intermediate.")
    
    ##more modifications, filter by tissue type
    parser.add_option('--expressedproteins',type='string',default='',dest='expr_prots',help='OPTIONAL: filters protein interaction network by protein identifiers provided in list. Best if drawn from mRNA Expression data for specific tissue, but also can be grabbed from gene atlas list of tissue-specific proteins')

    parser.add_option('--debug',action='store_true',default=False, dest='debug',help='Set to get extra information printed to files')
    ####################################################################

#    parser.add_option('--de_file', dest='de_file', type='string', default='', help='OPTINAL: tab delimited file containing differential expression for all types of genes (whether mRNAs,or other proteins). Could be the output of a genome-wide mRNA expression experiment. First column is gene name, and 2nd column is fold change in expression (or log fold change, or any other measure of the change in expression. Gene names should be same as in the output file, before identifier matching.')

#    parser.add_option('--de_cap',dest='de_cap',type='string',default='sink',help='OPTIONAL: determines how to use differential expression as capacities.  DEFAULT setting is \'sink\' to use as only capacities from mRNA to sink.  When set to \'source\' SAMNet will use differential expression capacities on source to genetic hit data, and if set to \'all\' SAMNet will use differential expression capacities on all nodes')

    parser.add_option('--amplPath',dest='ampl_path',type='string',default='/net/dorsal/apps/ampl/',help='Path to local ampl executable')
    (options, args) = parser.parse_args()

    print 'Running SAMNet...'

    ##format an options string to broadcast specifics as we process the input

    #######################PHENOTYPIC/INDIRECT/PROTEIN and DIRECT/miRNA weights#######################3
    print '.............upstream inputs..................'
    input_type='indirect (i.e. genetic)'
    if options.direct_weights!='':
        input_type=input_type+' and direct'
        direct_weights=parseIn.get_direct_target_weights(parseIn.multiple_args_into_one_dict(options.direct_weights,options.treat_names))
    else:
        direct_weights={}

    ##now process indirect weights - those from source
    if options.prot_weights!='':
        print 'This option is deprecated, please use --proteinWeightsByComm'
        indirect_weights=parseIn.get_weights_phen_source(parseIn.multiple_args_into_one_dict(options.prot_weights,options.treat_names))
        tn=re.split(',',options.treat_names)
   #     if len(tn)>0:
   #         tnames=[b+'_treatment' for b in tn]

    elif options.prot_weights_comm!='': 
        # using one comprehensive phos files for all commodities
        lf,tn,expprots=parseIn.by_comm_into_one_dict(options.prot_weights_comm,doUpper=(options.updateIds in ['human','mouse','yeast']))
        indirect_weights=parseIn.get_weights_phen_source(lf)
 
    else:
        sys.exit("Need some protein weight inputs!")

    ##now process treatment weights -- either costs in single commodity case, or amount of flow in multi
    treat_weights={}
    if len(options.treat_weights)>0:
        weight_vals=re.split(',',options.treat_weights)
    
#        if len(tn)==0:
#            tn=['treatment '+i for i in range(weight_vals)]
        for w in range(len(tn)):
            treat_weights[tnames[w]]=float(weight_vals[w])

    input_str='Processing '+input_type+' weights on interaction network'
    print input_str
    
    ################################PPI NETWORK######################################################
    print '...............protein-protein interaction network...................'
    ppi_str='With a gamma of '+str(options.gamma)+' and cap of '+str(options.cap)

    #PPI comes as a pkl graph with weights already incorporated in the edge information, so no edge addition necessary
    #remember the PPI graph (has information from multiple graphs given from user input)
    PPI_with_weights = parseIn.get_ppi_network(options.PPIfile,doUpper=(options.updateIds in ['human','mouse','yeast']))#networkx.read_gpickle(options.PPIfile)
    PPI_with_weights = networkx.DiGraph(PPI_with_weights)
#    if options.updateIds!='':
    ppi_str=ppi_str+'.  Using '+options.updateIds+' identifiers'
 
                
    if options.expr_prots!='':
        ppi_str=ppi_str+' and filtering by '+os.path.basename(options.expr_prots)
        PPI_with_weights=parseIn.filter_ppi(PPI_with_weights,open(options.expr_prots,'r'))

    expr_prots=PPI_with_weights.nodes()

 #   PPI_with_weights = multiple_args_into_one_list(options.PPIfile,True)
    print ppi_str

        
    ############################TF->mRNA network###################################################
    print ".............protein-DNA and transcriptional data..................."

    if options.tf_mrna_weights!='' and not options.noTfs:
        tf_string='Using tf-mRNA interactions from '+os.path.basename(options.tf_mrna_weights)
        addm=True
    else:
        tf_string='Not using any tf-MRNA interactions'
        addm=False

    if options.trafile!='':
        print "This option is depracated, plase use --traByComm"
        tradata = parseIn.multiple_args_into_one_dict(options.trafile,options.treat_names)
        weights_mRNA_to_sink = tfNetwork.get_weights_mRNA_sink(tradata,'foldchange',options.upordown,addMrna=addm)
        print 'mRNA data from: '+','.join(weights_mRNA_to_sink.keys())

    elif options.trafile_comm!='':
        tradata,ttn,expmrna = parseIn.by_comm_into_one_dict(options.trafile_comm,doUpper=(options.updateIds in ['human','mouse','yeast']))
        weights_mRNA_to_sink = tfNetwork.get_weights_mRNA_sink(tradata,'foldchange',options.upordown,addMrna=addm,doUpper=(options.updateIds in ['human','mouse','yeast']))
        print 'mRNA data from: '+','.join(weights_mRNA_to_sink.keys())

    elif options.raw_tf_data!='':
        print "This option is depracated, just use --noTfs flag!"
        tradata = parseIn.multiple_args_into_one_dict(options.raw_tf_data,options.treat_names)
        weights_mRNA_to_sink = tfNetwork.get_weights_mRNA_sink(tradata,'foldchange',options.upordown,addMrna=addm)
        print 'TF node data from: '+','.join(weights_mRNA_to_sink.keys())


    #remember the transcriptional graph
    graph_tr = networkx.DiGraph()
    if not options.noTfs:
        graph_tr = tfNetwork.make_tf_data_into_network(options.tf_mrna_weights,addmrna=addm,expressed_prot_list=expr_prots+expmrna,doUpper=(options.updateIds in ['human','mouse','yeast']))

    
    print tf_string

        
    #NAME OF OUTPUT FILE

    #directory
    Dirname = os.path.dirname(options.outputfile)    

    #name of file without directory
    #gene = '/'+os.path.basename(options.outputfile)+'_gam_'+options.gamma+options.foldtra
    gene=os.path.basename(options.outputfile)


    diff_ex_dict={}
#   if opts.de_file!='':
#        for row in open(opts.de_file,'r').readlines():
#            arr=row.strip().split('\t')
#            diff_ex_dict[arr[0]]=float(arr[1])
#    print final_weights
    ############################Hierarchical Network capacities######### ###########################
    node_caps=defaultdict(dict)## dictionary of all nodes for each commodity
    
    
    run_rn(PPI_with_weights,indirect_weights,direct_weights,graph_tr,mrna_weights=weights_mRNA_to_sink,output=options.outputfile,updateIds=options.updateIds,cap=options.cap,gamma=options.gamma,solver=options.solver,ampl_path=options.ampl_path,debug=options.debug,noTfs=options.noTfs,treat_weights=treat_weights,diff_ex_vals=diff_ex_dict,de_cap='sink',doMCF=len(indirect_weights)>1,makeCombined=options.make_combined,node_caps=node_caps,add_in_hier=options.hier_cap)

def run_rn(PPI_with_weights,indirect_weights,direct_weights,graph_tr,mrna_weights,output,updateIds='',cap=0.7,gamma=8,solver='loqo',ampl_path='',debug=False,noTfs=False,treat_weights=dict(), diff_ex_vals=dict(),de_cap='sink',doMCF=False,makeCombined=False,node_caps={},add_in_hier=True,sinkGamma=False):
#    print de_cap
    
    '''
    Ridiculously annoying function designed so SAMNet can also be run externally with pre-processed data supplied
    PPI_with_weights: networkx digraph object with string identifiers, already pre-filtered for expressed proteins
    indirect_weights**: dictionary of dictionaries of protein indentifiers as keys with values being their weights
    direct_weights:** dictionary of dictionaries (empty if only indirects are used) containing gene names as keys with values as their weights
    graph_tf: networkx digraph of transcriptional network with weights on the edges
    mrna_weights**: dictionary of mRNA weights to sink, or tf weights if noTfs flag was used
    output: outputstring to be used for files
    updateIds: species with which to update the protein network
#    expr: name of expression file for post-processing
    cap: value at which to cap protein interactions
    gamma: default gamma
    solver: default solver
    debug: set to True for extra info

    noTfs: set to true if the transcriptional network should be ignored and mrnaweights be linked directly to ppi
    diff_ex_vals: dictionary of differential expression values
    de_cap: set to 'source' or 'all' if you want to add diff_ex_vals-based capacities
    add_in_hier: Set this flag to set hierarchical capacities (1/10^x where x is unweighted sp to source), will populate node_caps for you
    ** -> These dictionaries are modified by the method, so be sure to pass in copies.
nn    
    Note on identifiers: It is important that all identifiers match.  SAMNet assumes that mRNA nodes are a unique set from the protein interactions.  In practice, protein identifiers (including indirect targets of miRNA) are in STRING identifier format and mRNA are in gene name with 'mrna' appended to the end. 
    '''
    
    #if not doMCF:
    #    doMCF=len(indirect_weights)>1 ##just double check!!
    #doMCF=True
   # print indirect_weights.keys()
#    if doMCF and len(node_caps)>0:
#        print 'Selected to run in multi-commodity mode with hierarchical capacities.  This has not been tested yet!'
    #EXPERIMENTAL DATA TO INCLUDE IN THE ANALYSIS
    #the phenotypic data to be included into interactome (only the proteins that are already in interactome)
    #phenres is portion of phenotypic data that can be analyzed with the given interactome
    #for each protein in the phenotypic data, keep it in phenres if it is contained in the PPI
    ##NEW: or if its in the TF->mrna network!
    #phenres = [x for x in weights_source_to_phen.keys() if PPI_with_weights.has_node(x)]
    
    allcoms=[treat for treat in indirect_weights.keys() if treat in mrna_weights.keys()]
    print 'Using '+str(len(allcoms))+' commodities that are found in both '+str(len(indirect_weights))+' protein weights and '+str(len(mrna_weights))+' sink weights'
    
    ubcids=['icrogid:5822231', 'icrogid:1938627', 'icrogid:2399853', 'icrogid:2068280','UBC_HUMAN','UBC']
    innet=[u for u in ubcids if u in PPI_with_weights.nodes()]
    print 'Removing '+str(len(innet))+' UBC identifiers from network'
    for i in innet:
        PPI_with_weights.remove_node(i)

    ##now load up the identifier matching
    if updateIds.lower()=='mouseiref':
        idfname=os.path.join(id_directory,'mouse_genename_to_9606mitabiref.pkl')
    else:
        idfname=os.path.join(id_directory,'humanUniprotHugoEntryMapping.pkl')
    if updateIds.lower()=='none':
	idfname=''
	idmap={}
    else:
        idmap=pickle.load(open(idfname,'r'))

    phenres = []
    idswap=defaultdict(list)## keep track of any ids we swapped so we can inform user about failed mappings
    if(len(indirect_weights)>0):
        for treat in indirect_weights.keys():
            if treat not in allcoms:
                indirect_weights.pop(treat,None)
                continue
            #print '--phen---'+treat
            for x in indirect_weights[treat].keys():#change for multisource
                if PPI_with_weights.has_node(x) or graph_tr.has_node(x):##forgot to check graphTr!
                    if x not in phenres:
                        phenres.append(x)
                elif x in idmap.keys():
                    for n in idmap[x]:
                        if PPI_with_weights.has_node(n) or graph_tr.has_node(n):
                            #add new node to list
                            indirect_weights[treat][n]=indirect_weights[treat][x]
                            if n not in phenres:
                                phenres.append(n)
                            idswap[n].append(x)##map replaced value with original
                            indirect_weights[treat][x]=0.0#zero this one out!
                else:
                    indirect_weights[treat][x]=0.0##zero this weight so it's not tallied..
                    if len(diff_ex_vals)>0 and treat in diff_ex_vals.keys():
                        diff_ex_vals[treat][x]=0.0
    print "how many proteinWeights are in interactome:"
    print len(phenres)

    count_phen=0
    count_indirect=0

    phensInInteractome=open(output+'proteinWeightsInInteractome.txt','w')
    phensInInteractome.write('Commodity\tFoundInNetwork\tIdentifier\n')
    for treat in indirect_weights.keys():
        for x in indirect_weights[treat].keys():
            if x in idswap.keys():
                val=x+' ('+','.join(idswap[x])+')'
            else:
                val=x
            if x not in phenres:
                count_indirect+=1
                phensInInteractome.write(treat+'\tNotInInteractome\t'+val+'\n')
            else:
                count_phen+=1
                phensInInteractome.write(treat+'\tInInteractome\t'+val+'\n')
    phensInInteractome.close()

    
    if (debug):
        network_inclusion_details=open(output+'input_included_in_final_network.txt','w')
        network_inclusion_details.write(str(count_phen)+' Protein Inputs present in interactome \n')
        network_inclusion_details.write(str(count_indirect)+' Indirect inputs present in interactome \n')
    
    directres = []
    if len(direct_weights)>0 and rawTf=='':
        for treat in direct_weights.keys():##changed for multi source
            if treat not in allcoms:
                direct_weights.pop(treat,None)
                next
            #print '--mirna---'+treat
            for x in direct_weights[treat].keys():
                if graph_tr.has_node(x):
                    if x not in directres:
                        directres.append(x)
                else:
                    direct_weights[treat][x]=0.0##zero it out so it's not totalled

        print "how many direct (e.g. miRNA) weights are in network"
        print len(directres)

    ##the transcriptional data to be included in the graph -
    ##only the mRNA which is contained in the transcriptional network
    trares = []
    if noTfs:  ##theis just means that mRNA weights are actually tf weights
        for treat in mrna_weights.keys():
            if treat not in allcoms:
                mrna_weights.pop(treat,None)
                next
            for x in mrna_weights[treat].keys():
                if PPI_with_weights.has_node(x):
                    if x not in trares:
                        trares.append(x)

                ##commented this on 2012-12-06 -- 
                elif PPI_with_weights.has_node(re.sub('mrna','',x)):
                    PPI_with_weights.add_edge(re.sub('mrna','',x),x,weight=0.00001)
                    if x not in trares:
                        trares.append(x)
                else:
                    mrna_weights[treat][x]==0.0
        print "how many direct sink weights are linked to the protein interaction network"
        print(len(trares))
    else:
        for treat in mrna_weights.keys():
            if treat not in allcoms:
                mrna_weights.pop(treat,None)
                next
            #print '--mrna--'+treat
            for x in mrna_weights[treat].keys():
                if graph_tr.has_node(x):
                    if x not in trares:
                        trares.append(x)
                else:
                    mrna_weights[treat][x]=0.0
        print "how many genes are linked to the TFs in the transcriptional network"
        print len(trares)

    count_tra=0
    trasInInteractome=open(output+'ExprInInteractome.txt','w')    
    trasInInteractome.write('Commodity\tFoundInNetwork\tIdentifier\n')

    for treat in mrna_weights.keys():
        for x in mrna_weights[treat].keys():
            if x in trares:
                trasInInteractome.write(treat+'\tInInteractome\t'+re.sub('mrna','',x)+'\n')
                count_tra+=1
            else:
                trasInInteractome.write(treat+'\tNotInInteractome\t'+re.sub('mrna','',x)+'\n')

    trasInInteractome.close()

    if(debug):
        network_inclusion_details.write(str(count_tra)+' Transcriptional mRNA nodes linked to TFs in interactome \n')

    #ADD SOURCE AND SINK TO PPI
    
    #include S1 and T1 (renamed from S and T to avoid conflicts),
    #the transcriptional and the genetic data into the interactome
    #this operation is possible only if the transcriptional and the genetic sets are not empty
    if (len(trares)>0) and (len(phenres)>0):
        #connect pehnotypic data (proteins) to source and add the corresponding weights
        #treat_weights is used differently in single vs multicommodity implementations
        if(len(treat_weights)==0):
            treat_weights=dict()
            for treat in indirect_weights.keys():
                print treat
                treat_weights[treat]=1.0
        #else:
        #    print treat_weights

        source='S1'
        sink='T1'

        
        #connect transcriptional data (mRNA) to the sink and add the corresponding weights
        for treat in mrna_weights.keys():
            PPI_with_weights.add_edge(treat+'_sink',sink,{'weight':treat_weights[treat]})
            for mrna_node in mrna_weights[treat].keys():
                weight = mrna_weights[treat][mrna_node]
              #  print treat, mrna_node,str(weight)
                if weight!=0.0:##this means that it is connected to the TF network
                    PPI_with_weights.add_edge(mrna_node,treat+'_sink',{'weight':weight})

        if len(directres)>0: ##if we are using miRNA data:
        ##ADD treatment node!!
            for treat in direct_weights.keys():
                PPI_with_weights.add_edge(source,treat,{'weight':treat_weights[treat]})
            ##then addd the direct to mrna
                for direct_node in directres:
                    if direct_node in direct_weights[treat].keys(): #double check now
                        weight=direct_weights[treat][direct_node]#/sum(direct_weights[treat].values())
                        if weight!=0.0:
                          #  print treat,direct_node,str(weight)
                            PPI_with_weights.add_edge(treat,direct_node,{'weight':weight})
        #now it will connect phens to indirect if phen_source was redone
        for treat in indirect_weights.keys():
            PPI_with_weights.add_edge(source,treat,{'weight':treat_weights[treat]})
            for prot in phenres:
                if prot in indirect_weights[treat].keys():##add double check
                    weight = indirect_weights[treat][prot]#/sum(indirect_weights[treat].values())
                    if weight!=0.0:
                    #    print treat,prot,str(weight)
                        PPI_with_weights.add_edge(treat,prot,{'weight':weight})
        


#        if not doMCF:# and len(directres)==0:
            ##this means that we're running the original RN
#            source=indirect_weights.keys()[0]
#        if len(mrna_weights)==1:
#            sink=mrna_weights.keys()[0]+'_sink'


        ##wait to add tf-> mrna weights until we know which tfs are in the network
     #ADD TRANSCRIPTIONAL DATA (MRNA NODES) TO PPI
    
    #transcription factors are already in the interactome
    #only the mRNA from the transcriptional data is inserted into the PPI, the rest of the mRNA (which would also be connected to the existing transcriptional factors is not relevant for the solution (so we don't add it, so we don't waste memory)
        if not noTfs:
            for mrna_node in trares:
                #find all nodes (TF) that point to that mRNA in the transcriptional network
                neigh=graph_tr.predecessors(mrna_node)
                #select neighbours of that mRNA that are not themselves mRNAs (so, they are TF involved in the production of the mRNA)
                neigh1=[x for x in neigh if 'mrna' not in x]
                #connect each of these neighbours to the mRNA and add the corresponding
                #weight from the transcriptional networ
                for tf in neigh1:
                    if PPI_with_weights.has_node(tf):
                        weight = graph_tr.get_edge_data(tf,mrna_node)['weight']
                        if weight>0.0:
                            PPI_with_weights.add_edge(tf, mrna_node,{'weight':weight})
                         

        ##only add capacities on targets if 
        if len(mrna_weights)>1:# and not doMCF: ##parameter not used for MCF
            target_cap=True
        else:
            target_cap=False
        
        ##now add in hierarchical capacities if we didn't already do so
        if add_in_hier and len(node_caps)==0:
            print '...............Setting capacities by distance from treatment...............'
            ##create shortest-distance network for each commodity to avoid unecessary shortcuts
            for node in PPI_with_weights.nodes():
                for treat in indirect_weights.keys():
                    if node!=treat and node!=treat+'_sink' and node!=sink:
                        spdist=1.0
                        try:
                            spdist=networkx.shortest_path_length(PPI_with_weights,treat,node,None)
                        except networkx.exception.NetworkXNoPath:
                            continue
                        node_cap=power(10.0,float(spdist)*-1.0)
                        node_caps[treat][node]=node_cap

        single_comms=[]
        sources,sinks=[],[]
        doMCF=True
        if doMCF:
            ##first write single commodity files
            if makeCombined:
                print 'Writing single commodity version of SAMNet files'
                for treat in indirect_weights.keys():
                    print 'Writing '+treat+' files'
                    if '/' in output: ##added to commodate non-directory paths
                        tmp_output=os.path.dirname(output)+'/'+treat+'_ONLY_'+os.path.basename(output)
                    else:
                        tmp_output=treat+'_ONLY_'+os.path.basename(output)
                        
                    wf.write_all_files(PPI_with_weights,trares,phenres,directres,tmp_output,treat,treat+'_sink',cap,str(gamma),solver,usetargetcapacity=target_cap,diff_ex_vals=diff_ex_vals,de_cap=de_cap,node_caps=node_caps,debug=debug,sinkGamma=sinkGamma)##DEFAULT is to add target capacities if we're not using direct/indirect responses, this should be changed
                    sources.append(treat)
                    sinks.append(treat+'_sink')
                    single_comms.append(tmp_output)
            ##now run multi commodity
            print "Writing MCF version of SAMNet files"
            output+='multiComm'
#            wf.write_mcf_files(PPI_with_weights,trares,phenres,directres,output,source,sink,cap,gamma,solver,usetargetcapacity=target_cap,diff_ex_vals=diff_ex_vals,de_cap=de_cap,node_caps=node_caps,debug=debug)##DEFAULT is to add target capacities if we're not using direct/indirect responses, this should be changed
            wf.write_mcf_files(PPI_with_weights,trares,phenres,output,source,sink,cap,str(gamma),solver,usetargetcapacity=target_cap,diff_ex_vals=diff_ex_vals,de_cap=de_cap,node_caps=node_caps,debug=debug,sinkGamma=sinkGamma)##DEFAULT is to add target capacities if we're not using direct/indirect responses, this should be changed
      #  else:
#            wf.write_all_files(PPI_with_weights,trares,phenres,directres,output,source,sink,cap,gamma,solver,usetargetcapacity=target_cap,diff_ex_vals=diff_ex_vals,de_cap=de_cap,node_caps=node_caps,debug=debug)##DEFAULT is to add target capacities if we're not using direct/indirect responses, this should be changed

       #     wf.write_all_files(PPI_with_weights,trares,phenres,output,source,sink,cap,str(gamma),solver,usetargetcapacity=target_cap,diff_ex_vals=diff_ex_vals,de_cap=de_cap,node_caps=node_caps,debug=debug,sinkGamma=sinkGamma)##DEFAULT is to add target capacities if we're not using direct/indirect responses, this should be changed
        
        #execute loqo
        single_comms.append(output)
        sources.append(source)
        sinks.append(sink)
        for i in range(len(single_comms)):
            out=single_comms[i]
            source=sources[i]
            sink=sinks[i]
            cmd=os.path.join(ampl_path,'ampl')
            args=out+".ampl"
    #cmd='ampl '+ Dirname+gene+".ampl"
            print 'Running '+solver+': '+time.asctime(time.localtime())
            retcode = subprocess.call(cmd+' '+args,shell=True)
            if retcode!=0:
                nr=subprocess.call('ampl_lic stop',shell=True)
                ##sometimes we have trouble getting license, try again
                retcode = subprocess.call(cmd+' '+args,shell=True)

            #            os.system(cmd)
            print 'Finished '+solver+': '+time.asctime(time.localtime())+' with return code '+str(retcode)
        ##call the post processing module to update the files
            if os.path.exists(out+'.txt'):
                if out==output: ##then we respect the original assignment
                    flow,node_flow,comm_flow,phens,prots,tfs,mrnas=post.process_output(out,source,sink, idfname,debug,diff_ex_vals,doMCF)
                else: ##otherwise we need to treat this as a non mcf
                    flow,node_flow,comm_flow,phens,prots,tfs,mrnas=post.process_output(out,source,sink, idfname,debug,diff_ex_vals,False)
                if flow>0.0:
            ##now output expression analysis
                    basename=out
                    noa='.noa'
                    if idfname!='' and idfname.lower()!='none':
                        noa='_symbol.noa'
            else:
                print 'Flow was 0, no output file created' ##NEED ERROR PROCESSING HERE

                


        if not debug:
            print 'Now removing .dat and .txt files to save space'
            for out in single_comms:
                try:
                    retcode=subprocess.call('rm '+out+'.txt '+out+'.dat',shell=True)
                except IOerror:
                    retcode=subprocess.call('rm '+out+'.dat',shell=True)
        if makeCombined and doMCF:
            single_comms=single_comms[0:len(single_comms)-1]
            edas=[f+'_ppi_attributes_symbol.eda' for f in single_comms]
#            if updateIds!='':
#                combine_single_flows_to_make_multi(edas,orig_output=output)
            edas=[f+'_ppi_attributes.eda' for f in single_comms]
            combine_single_flows_to_make_multi(edas,orig_output=output)
            for to in single_comms:
                retcode=subprocess.call('rm '+to+'*')

        return flow,phens,prots,tfs,mrnas
    #if either the transcriptional or the genetic set are empty, the program doesn't do anything
    else:
        print ""
        print "empty phenotypic or transcriptional response"



def combine_single_flows_to_make_multi(filename_list,orig_output,collapse_edges=False,ismcf=True):
    '''
    takesa  list of edge attribute files and combines them into sif and eda file...
    copied from ../chemoExpr/bin/compare_mcf_with_single.py
    Then modified to include a second sif file, this time with the results of a combined MCF/merged network
    Last argument collapses edges into a single edge, weighted by the fraction of commodities selecting that edge
    '''
    #create all 5 files
    final_mcf_sif=[]
    final_mcf_edge_comm=['EdgeCommodity\n']
    final_mcf_edge_type=['Interaction Type\n']
    
    final_mcf_node_type=['NodeType\n']
    final_mcf_node_flow=['Node Flow\n']
    final_mcf_node_comm_flow=[]


    ##dictionaries to handle new flow values
    node_flow_dict={}
    ##need to create dictionary to hold all info
    comm_flow_dict=defaultdict(dict)##indexed by gene, then commodity
    
    ##also open original 5 files for combined
    if 'symbol' in filename_list[0]:
        sym='_symbol'
    else:
        sym=''
#    if ismcf:
#        mc='_mcfs'
#    else:
#        mc='_all'

    sif_output=orig_output+'_mcfs'+sym+'.sif'
    edge_type_output=orig_output+'_edge_type'+sym+'.eda'
    edge_comm_output=orig_output+'_edge_commodity'+sym+'.eda'
    
    node_type_output=orig_output+'_node_type'+sym+'.noa'
    node_flow_output=orig_output+'_node_flow'+sym+'.noa'
    node_comm_output=orig_output+'_node_comm_flow'+sym+'.noa'


    ##dictionaries to handle existing flow values or counts if collapse_edges=T
    combined_node_flow_dict={}
    combined_comm_flow_dict=defaultdict(dict)

    ###check if files exist, then open them here
    if os.path.exists(sif_output):
        combined_mcf_sif=open(sif_output,'r').readlines()
    else:
        combined_mcf_sif=[]
        print 'no file: '+sif_output
        
    if os.path.exists(edge_comm_output):
        combined_mcf_edge_comm=open(edge_comm_output,'r').readlines()
    else:
        combined_mcf_edge_comm=['EdgeCommodity\n']
        print 'no file: '+edge_comm_output
        
    if os.path.exists(edge_type_output):
        combined_mcf_edge_type=open(edge_type_output,'r').readlines()
    else:
        combined_mcf_edge_type=['Interaction Type\n']
    
        
    if os.path.exists(node_type_output):
        combined_mcf_node_type=open(node_type_output,'r').readlines()
    else:
        combined_mcf_node_type=['NodeType\n']
        print "no file: "+node_type_output
        
    if os.path.exists(node_flow_output):
        combined_mcf_node_flow=open(node_flow_output,'r').readlines()
        for row in combined_mcf_node_flow[1:]:
            gene,val=row.strip().split(' = ')
            combined_node_flow_dict[gene]=float(val)
    else:
        print "no file: "+node_flow_output
            

    ##now re-assign        
    combined_mcf_node_flow=['Node Flow\n']

    ##get original flow data for each node
    if os.path.exists(node_comm_output) and ismcf:
        combined_mcf_node_comm_flow=open(node_comm_output,'r').readlines()
        comms=[s.strip() for s in combined_mcf_node_comm_flow[0].split('\t')[1:]]
        print comms
        for row in combined_mcf_node_comm_flow[1:]:
            arr=row.strip().split('\t')
            commvals={}
            for c in range(len(comms)):
                combined_comm_flow_dict[arr[0]][comms[c]]=arr[c+1]
#    else:
#        combined_mcf_node_comm_flow=[]
    
    #list of commodities processed, use this for counts if collapsing
    comlist=[]

    #create edge/flow dictionaries instead of writing to file
    sif_file_list=[]
    edge_comm_file_dict={}
    edge_type_file_dict={}
    
    for f in filename_list:
        ##break up filename
        if not os.path.exists(f):
            print 'File does not exist, likely no flow for this network:\n'+f
            continue
        rows=open(f,'r').readlines()[1:]
        f=os.path.basename(f)
        if 'Yeast' in f:
            commname=re.sub('Yeast','',f.split('_RN_treatment_ONLY')[0])
        else:
            commname=f.split('REMOVED')[0].split('_')[-2]
        comlist.append(commname)
        print 'Processing commodity '+commname
        flowvals=recalc_node_flow(rows,collapse_edges)
#        flowlist=[]
        typelist=[]        
        for row in rows[1:]:
            #print row
            p1,i_type,p2,eq,flow=row.strip().split()
            i_type=re.sub('\(','',re.sub('\)','',i_type))
            if i_type not in typelist:
                typelist.append(i_type)

            if collapse_edges and ismcf:
                commname=i_type
            #add edge to list
            edge=p1+'\t'+commname+'\t'+p2+'\n'
            if edge not in sif_file_list:
               # print edge
                sif_file_list.append(edge)

            ##now add counts attributes
            edge_a=p1+' ('+commname+') '+p2
            if edge_a not in edge_comm_file_dict.keys():
                edge_comm_file_dict[edge_a]=1.0
            else:
                edge_comm_file_dict[edge_a]=edge_comm_file_dict[edge_a]+1.0

            if edge_a not in edge_type_file_dict.keys():
                if collapse_edges:
                    edge_type_file_dict[edge_a]=flow.strip()
                else:
                    edge_type_file_dict[edge_a]=re.sub('\(','',re.sub('\)','',i_type))
            
            #            if p2 not in flowlist:
            #                flowlist.append(p2)
            if p2 in node_flow_dict.keys():
                node_flow_dict[p2]+=flowvals[p2]
            else:
                node_flow_dict[p2]=flowvals[p2]

            if p2 in combined_node_flow_dict.keys():
                combined_node_flow_dict[p2]+=flowvals[p2]
            else:
                combined_node_flow_dict[p2]=flowvals[p2]

            ##now do the commodity flow
            if p2 in comm_flow_dict.keys() and commname in comm_flow_dict[p2].keys():
                comm_flow_dict[p2][commname]+=flowvals[p2]
            else:
                comm_flow_dict[p2][commname]=flowvals[p2]
            if p2 in combined_comm_flow_dict.keys() and commname+'_altered' in combined_comm_flow_dict[p2].keys():
                combined_comm_flow_dict[p2][commname+'_altered']+=flowvals[p2]
            else:
                combined_comm_flow_dict[p2][commname+'_altered']=flowvals[p2]

#            final_mcf_node_flow.append(p2+' = '+str(flowvals[p2])+'\n')
            node_type_list=[]
            if p2 not in node_type_list:
                if p1 in ['S1','arsenic','copper','cadmium','chromium','mercury','silver','zinc','Zeb1','Snail','Tgfb','Fixed']:
                    node_type_list.append(p2)
                    final_mcf_node_type.append(p2+' = phenotypic\n')
                    combined_mcf_node_type.append(p2+' = phenotypic\n')
                    
            if p1 not in node_type_list:
                if p2=='T1' or 'sink' in p2:
                    node_type_list.append(p1)
                    final_mcf_node_type.append(p1+' = mrna\n')
                    combined_mcf_node_type.append(p1+' = mrna\n')
                if 'H3K' in p1:
                    if 'H3K4me3' in p1:
                        histone='H3K4me3'
                    elif 'H3K27ac' in p1:
                        histone='H3K27ac'
                    elif 'H3K36me3' in p1:
                        histone='H3K36me3'
                    else:
                        histone=''
                    node_type_list.append(p1)
                    final_mcf_node_type.append(p1+' = '+histone+'transcriptionfactor\n')
                    combined_mcf_node_type.append(p1+' = '+histone+'transcriptionfactor\n')


    #first handle node flow commodities
    final_mcf_node_flow.extend([node+' = '+str(node_flow_dict[node])+'\n' for node in node_flow_dict.keys()])
    combined_mcf_node_flow.extend([node+' = '+str(combined_node_flow_dict[node])+'\n' for node in combined_node_flow_dict.keys()])

    ##then handle edge files
    ##add edge to sif files          
    final_mcf_sif=sif_file_list
    for edge in sif_file_list:
        vals=edge.split('\t')
        vals[1]=vals[1]+'_altered'
       #this seems wrong vals[1]=vals[2]+'_altered'
        combined_mcf_sif.append('\t'.join(vals))
    
    #add edge to attribute files
    for edge_a in edge_comm_file_dict.keys():
        cv=edge_comm_file_dict[edge_a]/len(comlist)
        et=edge_type_file_dict[edge_a]
        final_mcf_edge_comm.append(edge_a+' = '+str(cv)+'\n')
        final_mcf_edge_type.append(edge_a+' = '+et+'\n')
        
        edge_arr=edge_a.split()
        ##replaced arr with edge_arr below, will this break actual mcf code?  seems like a bug
        edge_arr[1]='('+re.sub('\(','',re.sub('\)','',edge_arr[1]))+'_altered)'
        combined_mcf_edge_comm.append(' '.join(edge_arr)+' = '+str(cv)+'\n')
        combined_mcf_edge_type.append(' '.join(edge_arr)+' = '+et+'\n')
    
    ##now handle node flow matrix file
    comms=comlist
    if collapse_edges and ismcf:
        comms=typelist
    final_mcf_node_comm_flow=['Node\t'+'\t'.join(comms)+'\n']
    for node in comm_flow_dict.keys():
        row=node
        for comm in comms:
            if comm in comm_flow_dict[node].keys():
                row+='\t'+str(comm_flow_dict[node][comm])
            else:
                row+='\t0'
        final_mcf_node_comm_flow.append(row+'\n')

    ##now make node flow commodity file for combined
    combined_comms=comms+[c+'_altered' for c in comms]
    print combined_comms
    combined_mcf_node_comm_flow=['Node\t'+'\t'.join(combined_comms)+'\n']
    for node in combined_comm_flow_dict.keys():
        row=node
        for com in combined_comms:
            if com in combined_comm_flow_dict[node].keys():
                row+='\t'+str(combined_comm_flow_dict[node][com])
            else:
                row+='\t0'
        combined_mcf_node_comm_flow.append(row+'\n')
        
    orig_output=os.path.dirname(orig_output)+'/MERGED_FROM_SINGLE'+re.sub('commFlow','',os.path.basename(orig_output))
    combined_output=os.path.dirname(orig_output)+'/MERGED_And_COMBINED'+re.sub('commFlow','',os.path.basename(orig_output))
    
    newfname=orig_output+'_mcfs'+sym+'.sif'
    newedaname=orig_output+'_edge_commodity'+sym+'.eda'
    new_type_edaname=orig_output+'_edge_type'+sym+'.eda'    
    newtypename=orig_output+'_node_type'+sym+'.noa'
    newflowname=orig_output+'_node_flow'+sym+'.noa'
    newnodeflowname=orig_output+'_node_comm_flow'+sym+'.noa'
  #  if not os.path.exists('combined_graphs'):
  #      os.system('mkdir combined_graphs')
    open(newfname,'w').writelines(final_mcf_sif)
    open(newedaname,'w').writelines(final_mcf_edge_comm)
    open(newtypename,'w').writelines(final_mcf_node_type)    
    open(newflowname,'w').writelines(final_mcf_node_flow)
    open(newnodeflowname,'w').writelines(final_mcf_node_comm_flow)
    open(new_type_edaname,'w').writelines(final_mcf_edge_type)

    ##now write combined files
    comb_newfname=combined_output+'_mcfs'+sym+'.sif'
    comb_newedaname=combined_output+'_edge_commodity'+sym+'.eda'
    comb_type_edaname=combined_output+'_edge_type'+sym+'.eda'
    comb_newtypename=combined_output+'_node_type'+sym+'.noa'
    comb_newflowname=combined_output+'_node_flow'+sym+'.noa'
    comb_newnodeflowname=combined_output+'_node_comm_flow'+sym+'.noa'
    
    open(comb_newfname,'w').writelines(combined_mcf_sif)
    open(comb_newedaname,'w').writelines(combined_mcf_edge_comm)
    open(comb_type_edaname,'w').writelines(combined_mcf_edge_type)
    open(comb_newtypename,'w').writelines(combined_mcf_node_type)    
    open(comb_newflowname,'w').writelines(combined_mcf_node_flow)
    open(comb_newnodeflowname,'w').writelines(combined_mcf_node_comm_flow)

def recalc_node_flow(eda_rows,not_flow=True):
    flowdict={}
    for row in eda_rows:
      #  print row
        p1,itype,p2,eq,flow=row.strip().split()
        if not_flow:
            flowdict[p2]=1.0
        elif p2 not in flowdict.keys():
            flowdict[p2]=float(flow)
        else:
            flowdict[p2]+=float(flow)
    return flowdict


if __name__=='__main__':
    main()
