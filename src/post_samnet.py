'''
post_samnet.py

SAMNet module designed to post-process the output of AMPL into cytoscape-readable files

Copyright (c) 2012 Sara JC Gosline
sgosline@mit.edu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE


'''
#import enrichmentStats as es
import idmatch_samnet as identifier_matching
import re,pickle,os
from collections import defaultdict

fpath=os.path.dirname(os.path.abspath( __file__ ))
#print fpath

id_directory=re.sub('src','lib',fpath)
#print id_directory


#write sif file
#removed chosen phens and chose tra file, one less to worry about
def write_sif_file(wholename, source, sink,node_flow,comm_flow,debug=False, de_file=None,mcf=False):
  #  print 'inside write sif'
    #print de_file
    file = open(wholename+'.txt', 'r')
    lines = file.readlines()
    siffile1 = open(wholename+'_no_mrna.sif', 'w')
    siffile2 = open(wholename+'_tf_2_mrna.sif', 'w')
    siffile3=open(wholename+'_all.sif','w')
    siffile4=open(wholename+'_mcfs.sif','w')

    ##create specific sif files for each commodity.  then can run in cytoscape with existing noa and eda files
    comm_sif_files=dict()
    comm_eda_files=dict()

    attrfile1=open(wholename+'_ppi_attributes.eda','w')
    attrfile1.write('FlowThroughEdge'+'\n')
    attrfile2=open(wholename+'_node_type.noa','w')
    attrfile2.write('NodeType'+'\n')
    attrfile3=open(wholename+'_node_flow.noa','w')
    attrfile3.write("Node Flow"+'\n')
    

    ##this is for MCF files only, its included as interaction type in responsenet
    attrfile5=open(wholename+'_edge_commodity.eda','w')
    attrfile6=open(wholename+'_edge_type.eda','w')
    attrfile7=open(wholename+'_node_comm_flow.noa','w')

    if(mcf):
        attrfile5.write('EdgeCommodity'+'\n')
        attrfile6.write('Interaction Type\n')
        attrfile7.write('Node\t'+'\t'.join(comm_flow.keys())+'\n')
        for n in node_flow.keys():
            nodestr=n
            for c in comm_flow.keys():
                if n in comm_flow[c].keys():
                    nc=comm_flow[c][n]
                else:
                    nc=0.0
                nodestr+='\t'+str(nc)
            attrfile7.write(nodestr+'\n')
        find=3#set index of flow in txt file
    else:
        find=2

    if len(de_file)>0:
        #add other attribute files for expression of genes
        attrfile4=open(wholename+'_DiffExpr.noa','w')
        attrfile4.write("DiffExpression"+"(class=Double)"+'\n')
    
    '''if (debug):
        print 'debug'
        phensInInteractome=open(wholename+'proteinWeightsInInteractome.txt','r')
        listProtWeights=[]
        #print phensInInteractome.readlines()
        a=phensInInteractome.readlines()
        print a
        for sth in a:
            
            listProtWeights.append(sth.strip('\n'))
        print len(listProtWeights)
        print listProtWeights
        protWeightsPicked=[]
    '''
    
    phens=set()
    other_prots=set()
    tfs=set()
    mrnas=set()
    dirmrnas=set()
    
    if(node_flow==-1):
        return phens,other_prots,tfs,mrnas

   # tfs,dirmrna,indirmrna,phens,flows=[],[],[],[],[]

#    if de_file!=None:
#        dExp_file=open(de_file,'r')
    for gene in de_file.keys():
#        gene=fields[0]
        de=de_file[gene]
#        de=float(fields[1])
        if de >=100:
            de=100
        elif de<=-100:
            de=-100
            #add differential expression
        attrfile4.write(gene+' = '+str(de)+'\n')

    ##keep set to properly annotate nodes that belong to multiple categories

    #indirmrna=set()
    flows=set()
    for line in lines:
        fields = line.split('\t')
        prot1 = fields[0].strip()
        prot2 = fields[1].strip()
        flow = str(fields[find].strip())
        if flow=='0.000000' or flow=='-0.000000':
            continue
        comm=''
        if mcf:
            comm=fields[2].strip()
        ##first check if prot1 is a source or multi source
        ##label this as pm interaction indicating non pp or pd

        #first evaluate phenotypic hits
        if (prot1==source or 'treatment' in prot1) and prot2!=sink:
            tprot1=re.sub('_treatment','',prot1)
            tprot2=re.sub('mrna','',prot2)
  #          tprot2=re.sub('_treatment','',prot2)
            siffile1.write(tprot1+'\t'+'pm'+'\t'+tprot2+'\n')
            siffile3.write(tprot1+'\t'+'pm'+'\t'+tprot2+'\n')
            siffile4.write(tprot1+'\t'+comm+'\t'+tprot2+'\n')
            attrfile1.write(tprot1+' (pm) '+tprot2+' = '+flow+'\n')
            attrfile5.write(tprot1+' ('+comm+') '+tprot2+' = '+flow+'\n')
            #if the phenotypic/mirna hit is indeed an mrna, make it a diamond
            if 'mrna' in prot2 and prot2 not in mrnas:
                mrnas.add(tprot2)
                #attrfile2.write(tprot2+' = '+'mrna'+'\n')
                #protein mrna
            attrfile6.write(tprot1+' ('+comm+') '+tprot2+' = pm\n')
            #otherwise, make it a square
#            elif tprot2 not in phens and tprot2 not in mrnas:
            phens.add(tprot2)

            if tprot2 not in flows and tprot2 in node_flow.keys():
                attrfile3.write(tprot2+' = '+str(node_flow[tprot2])+'\n')
                flows.add(tprot2)

        #if edge is within proteins, write protein interaction append flow info
        elif prot2!=sink and 'mrna' not in prot2 and 'mrna' not in prot1:
            siffile1.write(prot1+'\t'+'pp'+'\t'+prot2+'\n')
            siffile3.write(prot1+'\t'+'pp'+'\t'+prot2+'\n')
            siffile4.write(prot1+'\t'+comm+'\t'+prot2+'\n')
            attrfile1.write(prot1+' (pp) '+prot2+' = '+flow+'\n')
            attrfile5.write(prot1+' ('+comm+') '+prot2+' = '+flow+'\n')
            attrfile6.write(prot1+' ('+comm+') '+prot2+' = pp\n')
#            other_prots.add(prot1)
            other_prots.add(prot2)
            if prot2 not in flows and prot2 in node_flow.keys():
                attrfile3.write(prot2+' = '+str(node_flow[prot2])+'\n')
                flows.add(prot2)

        #if edge goes from protein (or direct) to mrna, label as mrna, add flow info
        elif 'mrna' in prot1:
            prot1=re.sub('mrna','',prot1)
            mrnas.add(prot1)
            prot2=re.sub('_treatment','',prot2)
            siffile2.write(prot1+'\t'+'pm'+'\t'+prot2+'\n')
            siffile3.write(prot1+'\t'+'pm'+'\t'+prot2+'\n')
            siffile4.write(prot1+'\t'+comm+'\t'+prot2+'\n')
            attrfile1.write(prot1+' (pm) '+prot2+' = '+flow+'\n')
            attrfile5.write(prot1+' ('+comm+') '+prot2+' = '+flow+'\n')
            attrfile6.write(prot1+' ('+comm+') '+prot2+' = pm\n')
        elif 'mrna' in prot2:
            prot2=re.sub('mrna','',prot2) ##this could induce loops

            mrnas.add(prot2)
            tfs.add(prot1)            
            siffile2.write(prot1+'\t'+'pd'+'\t'+prot2+'\n')
            siffile3.write(prot1+'\t'+'pd'+'\t'+prot2+'\n')
            siffile4.write(prot1+'\t'+comm+'\t'+prot2+'\n')
            attrfile6.write(prot1+' ('+comm+') '+prot2+' = pd\n')
            attrfile1.write(prot1+' (pd) '+prot2+' = '+flow+'\n')
            attrfile5.write(prot1+' ('+comm+') '+prot2+' = '+flow+'\n')
            if prot2 not in flows and prot2 in node_flow.keys():
                attrfile3.write(prot2+' = '+str(node_flow[prot2])+'\n')
                flows.add(prot2)
                

      #  elif prot2==sink and 'treatment' not in prot1::
      #      if 'mrna' not in prot1:# and prot1 not in tfs:
      #          tfs.add(prot1)
                #attrfile2.write(prot1+' = '+'transcriptionfactor'+'\n')

    #now go through all proteins in all sets and add attribute
    allprots=set()
    allprots.update(phens)
    allprots.update(mrnas)
    allprots.update(tfs)
    allprots.update(other_prots)
    for prot in allprots:
        attrstring=prot+' = '
        if prot in tfs:
            attrstring+='transcriptionfactor'
        if prot in mrnas:
            attrstring+='mrna'
        if prot in phens:
            attrstring+='phenotypic'
        attrfile2.write(attrstring+'\n')

    
    if (debug):
        print 'writing to details'
        network_output_details=open(wholename+'input_included_in_final_network2.txt','w')
        network_output_details.write(str(len(phens))+' Phenotypic included in final network')
        network_output_details.write(str(len(mrnas))+' Indirect mRNA included in final network')
#        network_output_details.write(str(len(dirmrna))+' Direct mRNA included in final network')
        network_output_details.write(str(len(tfs))+' TFs included in final network')

        '''
        network_output_details.write(str(len(protWeightsPicked))+' Total protWeights included in final network (also includes those not connected to source)')
        '''  
    return phens,other_prots,tfs,mrnas
            


##calculates incoming flow out of total flow for each node
def calculate_node_flow(lines,mcf):
    node_flow={}
    total=0.0
    comm_flow=defaultdict(dict)
    #TODO: calculate node flow individually for each commodity?
    if(mcf):
        find=3
    else:
        find=2
    
    for l in lines:
        l=re.split('\t',l.strip())
        ind=re.sub('mrna','',l[1])
        if float(l[find].strip())==0.0:
            continue
        if ind in node_flow.keys():
            node_flow[ind]=node_flow[ind]+float(l[find].strip())
        else:
            node_flow[ind]=float(l[find].strip())
        if mcf:
            comm=l[2].strip()
            if comm in comm_flow.keys():
                if ind in comm_flow[comm].keys():
                    comm_flow[comm][ind]+=float(l[find].strip())
                else:
                    comm_flow[comm][ind]=float(l[find].strip())
            else:
                comm_flow[comm][ind]=float(l[find].strip())
        total=total+float(l[find].strip())
    if(total==0):
        print "No flow, try increasing gamma"
        return {},{},0

    norm_flow={}
    for k in node_flow.keys():
        norm_flow[k]=node_flow[k]/total

    print 'Total flow:',total
    return norm_flow,comm_flow,total


def process_output(output_file,source='S', sink='T', species_name='',debug=False,de_file=None,mcf=False):
    '''
    Run the standard post-processing steps for responseNet
    '''
    
    if not os.path.exists(output_file+'.txt'):
        print 'Output file missing'
        return dict(),dict(),0.0,set(),set(),set()
   ##Calculate node flow for ranking of signaling proteins
    (node_flow,comm_flow,total)=calculate_node_flow(open(output_file+'.txt','r').readlines(),mcf)#returns a dictionary of node flow
    
        #calculate enrichment statistic if mRNA are used?

        #visualize
    if total==0.0:
        print 'No flow'
        return total,node_flow,comm_flow,set(),set(),set()
    phens,prots,tfs,mrnas=write_sif_file(output_file, source, sink,node_flow,comm_flow,debug,de_file,mcf)
            ##MODIFIED by SGOSLINE: added this to do identifer matching for the sif files
    
    if(species_name.lower==''):
        print 'No identifier matching, moving on...'
        idfile=''
    else:
        if(species_name.lower()=='mouse'):
            idfile=pickle.load(open(id_directory+'/10090protein.aliases.v9.0_geneName.pkl','r'))
        elif(species_name.lower()=='human'):
            idfile=pickle.load(open(id_directory+'/9606protein.aliases.v9.0_geneName.pkl','r'))
        elif(species_name.lower()=='yeast'):
            idfile=pickle.load(open(id_directory+'/4932protein.aliases.v9.0_geneName.pkl','r'))
        elif(species_name.lower()=='humaniref'):
            idfile=pickle.load(open(id_directory+'/9606mitab.01192011.uniq_miscore-localirefindex3-20110831.geneMapping.pkl','r'))
        elif(species_name.lower()=='mouseiref'):
            idfile=pickle.load(open(id_directory+'/mouse_genename_to_9606mitabiref.pkl','r'))
	else:
            idfile=''
    if idfile!='': 

        print "Matching identifiers"
        identifier_matching.parseSifFileFromStringToGeneName(open(output_file+'_all.sif','r'),output_file+'_all_symbol.sif',idfile)
        identifier_matching.parseSifFileFromStringToGeneName(open(output_file+'_mcfs.sif','r'),output_file+'_mcfs_symbol.sif',idfile)
        identifier_matching.parseSifFileFromStringToGeneName(open(output_file+'_no_mrna.sif','r'),output_file+'_no_mrna_symbol.sif',idfile)
        #also for the edge attribute files
        identifier_matching.parseAttrFileFromStringToGeneName(open(output_file+'_ppi_attributes.eda','r'),output_file+'_ppi_attributes_symbol.eda',idfile)        
        identifier_matching.parseAttrFileFromStringToGeneName(open(output_file+'_edge_commodity.eda','r'),output_file+'_edge_commodity_symbol.eda',idfile)
        identifier_matching.parseAttrFileFromStringToGeneName(open(output_file+'_edge_type.eda','r'),output_file+'_edge_type_symbol.eda',idfile)

        identifier_matching.pareTabFileFromStringToGeneName(open(output_file+'_node_comm_flow.noa','r'),output_file+'_node_comm_flow_symbol.noa',idfile)

##created new function for node attributes
        identifier_matching.parseNodeAttrFileFromStringToGeneName(open(output_file+'_node_type.noa','r'),output_file+'_node_type_symbol.noa',idfile,True)
        identifier_matching.parseNodeAttrFileFromStringToGeneName(open(output_file+'_node_flow.noa','r'),output_file+'_node_flow_symbol.noa',idfile,False)
        if len(de_file)>0:
            identifier_matching.parseNodeAttrFileFromStringToGeneName(open(output_file+'_DiffExpr.noa','r'),output_file+'_DiffExpr.noa',idfile,False)
    return total,node_flow,comm_flow,phens,prots,tfs,mrnas
