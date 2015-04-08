'''
writefiles_samnet.py
SAMNet module responsible for writing files to be used by AMPL


Copyright (c) 2012 Sara JC Gosline
sgosline@mit.edu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

'''

import collections,re,networkx

#write .mod file
def writechangeflow(wholename,gamma,sinkGamma=False):
    file=open(wholename+'changeflow.mod','w')
    file.write('set proteins;\n')
    file.write('set sourcesink;\n')
    file.write('set nodes = sourcesink union proteins;\n')
    file.write('set interactions within {nodes cross nodes};\n')
    file.write('set sink_interactions within {i in proteins,j in sourcesink};\n')
    file.write('set source_interactions within {i in sourcesink,j in proteins};\n')
    file.write('param cost {interactions} >=0;\n')
    file.write('param capacity {interactions} >=0;\n')
    file.write('var X {(i,j) in interactions}>=0, <=capacity[i,j];\n')
    file.write('minimize Total_Cost:')
    if sinkGamma:
        file.write('sum{(i,j) in interactions}-log(cost[i,j])*X[i,j]-sum{(i,j) in sink_interactions}'+gamma+'*X[i,j];\n')
    else:
        file.write('sum{(i,j) in interactions}-log(cost[i,j])*X[i,j]-sum{(i,j) in source_interactions}'+gamma+'*X[i,j];\n')
    file.write('subject to Kirkoff {k in proteins}: sum {(i,k) in interactions} X[i,k]=sum{(k,j) in interactions} X[k,j];\n')
    file.write('subject to sourcesinkcond : sum{(i,j) in source_interactions} X[i,j]=sum{(i,j) in sink_interactions} X[i,j];\n')
    file.close()

              


def write_mcf_changeflow(wholename,gamma,revGamma=False):
    file=open(wholename+'changeflow.mod','w')
    file.write('set proteins;\n')
    file.write('set source;\n')
    file.write('set sink;\n')
    file.write('set commodities;\n') ##added this
    file.write('set initnodes;\n')    ##need to add these
    file.write('set endnodes;\n') ##and these
    file.write('set termnodes  = initnodes union endnodes;\n') #created this
    file.write('set interiornodes = proteins diff termnodes;\n') #and this
    file.write('set sourcesink = source union sink;\n')
    file.write('set nodes = sourcesink union proteins;\n')
    file.write('set interactions within {proteins cross proteins};\n')
    file.write('set all_interactions within {nodes cross nodes};\n')

    file.write('set sink_interactions{commodities} within {i in proteins,j in sink};\n')
    file.write('set source_interactions{commodities} within {i in source,j in proteins};\n')

    file.write('param cost {all_interactions,commodities} >=0;\n') ##every interaction has a cost
    
    file.write('param capacity {all_interactions} >=0;\n') ##every interaction has a capacity, most are 1
    file.write('var X {all_interactions,commodities} >=0;\n')
    file.write('minimize Total_Cost:\n')
    ##if this flag is used, gamma constraint is applied to interactions with sink, not source!
    if revGamma:
        file.write('sum{k in commodities}sum{(i,j) in interactions}-log(cost[i,j,k])*X[i,j,k] + sum{k in commodities}sum{(i,j) in source_interactions[k]}-log(cost[i,j,k])*X[i,j,k] + sum{k in commodities}sum{(i,j) in sink_interactions[k]}-log(cost[i,j,k])*X[i,j,k] - sum{k in commodities}sum{(i,j) in sink_interactions[k]}'+gamma+'*X[i,j,k];\n')
    else:
        file.write('sum{k in commodities}sum{(i,j) in interactions}-log(cost[i,j,k])*X[i,j,k] + sum{k in commodities}sum{(i,j) in source_interactions[k]}-log(cost[i,j,k])*X[i,j,k] + sum{k in commodities}sum{(i,j) in sink_interactions[k]}-log(cost[i,j,k])*X[i,j,k] - sum{k in commodities}sum{(i,j) in source_interactions[k]}'+gamma+'*X[i,j,k];\n')
#        file.write('sum{k in commodities}sum{(i,j) in all_interactions}-log(cost[i,j,k])*X[i,j,k] - sum{k in commodities}sum{(i,j) in source_interactions[k]}'+gamma+'*X[i,j,k];\n')
#    file.write('sum{(i,j) in interactions,k in commodities}-log(cost[i,j,k])*X[i,j,k]-sum{(i,j) in source_interactions,k in commodities}'+gamma+'*X[i,j,k];\n')
#    file.write('subject to Kirkoff {k in proteins,c in commodities}: sum {(i,k) in all_interactions} X[i,k,c]=sum{(k,j) in all_interactions} X[k,j,c];\n')
    file.write('subject to Interior_Kirkoff {k in interiornodes,c in commodities}: sum {(i,k) in interactions} X[i,k,c]=sum{(k,j) in interactions} X[k,j,c];\n')
    file.write('subject to Source_Kirkoff {k in initnodes,c in commodities}: sum{(i,k) in source_interactions[c]} X[i,k,c]=sum{(k,j) in interactions} X[k,j,c];\n')
    file.write('subject to Sink_Kirkoff {k in endnodes,c in commodities}: sum{(i,k) in interactions} X[i,k,c]=sum{(k,j) in sink_interactions[c]} X[k,j,c];\n')
    file.write('subject to sourcesinkcond{k in commodities}: sum{(i,j) in source_interactions[k]} X[i,j,k]=sum{(i,j) in sink_interactions[k]} X[i,j,k];\n')
    file.write('subject to multi {(i,j) in all_interactions}: sum {k in commodities}X[i,j,k]<=capacity[i,j];\n')##add constraint here
    file.close()

              
#def write_mcf_datfile(big_PPI,trares,phenres,directres,outputfilename,source,sink,cap,usetargetcapacity,diff_ex_vals,de_cap,node_caps={},debug=False):
def write_mcf_datfile(big_PPI,trares,phenres,outputfilename,source,sink,cap,usetargetcapacity,diff_ex_vals,de_cap,node_caps={},debug=False,norm_weights_and_caps=False):
#    print 'Writing mcf file'
    '''
    This is similar to the writedat file with multiple sources and sinks, but 
    compresses all sources and sinks to individual commodities.  Costs from source to hits
    and from tra to sink vary based on commodity, but internal edge weights are the same
    '''

    
    print 'Writing network with '+str(big_PPI.number_of_nodes())+' nodes and '+str(big_PPI.number_of_edges())+' edges'
    for edge in big_PPI.edges():
        n1=edge[0]
        n2=edge[1]
        ew=float(big_PPI.get_edge_data(n1,n2)['weight'])
        if ew==0.0:
            big_PPI.remove_edge(n1,n2)

    print 'After zero-weight edges removed, has '+str(big_PPI.number_of_nodes())+' nodes and '+str(big_PPI.number_of_edges())+' edges'
    ##extract commodity weights:
    suff=''
    #suff=suff
    comm_weights={}
    for c in big_PPI.successors(source):
        comm_weights[re.sub(suff,'',c)]=big_PPI.get_edge_data(source,c)['weight']

   # print comm_weights

    directres=[]
    #first let's get a handle on the commodities
    com_sources=big_PPI.successors(source)
    commodity_names=[]

    for c in com_sources:
        commodity_names.append(re.sub(suff,'',c))
    dummy_nodes=commodity_names+[a+'_sink' for a in commodity_names]
    print "All commodities: "+','.join(commodity_names)

    #now collect dictionary of commodity weights
    commodity_source_weights=collections.defaultdict(dict)
    commodity_direct_weights=collections.defaultdict(dict)

    for c in commodity_names:
        corig=c.strip('\"')
        for neigh in big_PPI.successors(corig+suff):
            if neigh in phenres:
                commodity_source_weights[c][neigh]=comm_weights[c]*big_PPI.get_edge_data(corig+suff,neigh)['weight']
            if neigh in directres:
                commodity_direct_weights[c][neigh]=comm_weights[c]*big_PPI.get_edge_data(corig+suff,neigh)['weight']
            #    print commodity_source_weights


    commodity_sink_weights=collections.defaultdict(dict)
    for c in commodity_names:
        corig=c.strip('\"')
        for neigh in big_PPI.predecessors(corig+suff+'_sink'):
            commodity_sink_weights[c][neigh]=abs(big_PPI.get_edge_data(neigh,corig+suff+'_sink')['weight'])
#    print commodity_sink_weights

    file=open(outputfilename+'.dat','w')
    file.write('set proteins :=')
    '''
    include everything in proteins but the multi sources and multi sinks
    '''
    for p in big_PPI.nodes():
        #if p not in big_PPI.successors(source) and p not in big_PPI.predecessors(sink)
        if p not in commodity_names and p!=source and p!=sink:
            file.write(' \"'+p+'\"')
    file.write(';\n')

    file.write('set source := '+source+';\n')
    file.write('set sink := '+sink+';\n')
    file.write('set commodities := '+' '.join(commodity_names)+';\n')
    file.write('set initnodes :='+' '.join(['"'+a+'"' for a in phenres])+';\n')
    file.write('set endnodes :='+' '.join(['"'+a+'"' for a in trares])+';\n')
    file.write('set interactions :=')
    #get all edges here
    for two_edge_nodes in big_PPI.edges():
        n1=two_edge_nodes[0]
        n2=two_edge_nodes[1]
        ew=float(big_PPI.get_edge_data(n1,n2)['weight'])

        if n1!=source and n2!=sink and n1 not in dummy_nodes and n2 not in dummy_nodes and ew!=0.0:
            file.write('(\"'+n1+'\",\"'+n2+'\")') 
    file.write(";\n")


    
    file.write('set all_interactions :=')
    #get all edges here
    for two_edge_nodes in big_PPI.edges():
        n1=two_edge_nodes[0]
        n2=two_edge_nodes[1]
        ew=float(big_PPI.get_edge_data(n1,n2)['weight'])
        
        if n1!=source and n2!=sink and n1 not in dummy_nodes and n2 not in dummy_nodes and ew!=0.0:
            file.write('(\"'+n1+'\",\"'+n2+'\")') 
    for p in phenres: #assume this is the union of all phenotypic data
        file.write('('+source+',\"'+p+'\")')
    for d in directres: ##need to account for direct interactions
        file.write('('+source+',\"'+d+'\")')
    for t in trares:#assume this is the union of all expression data
        file.write('(\"'+t+'\",'+sink+')')
    file.write(';\n')


    for c in commodity_names:
        file.write('set source_interactions['+c+'] :=')
        corig=c.strip('\"')
        for neighbors in big_PPI.successors(corig+suff):
            if neighbors in phenres+directres:
                file.write('(\"'+source+'\",\"'+neighbors+'\")')
        file.write(';\n')

    for c in commodity_names:
        file.write('set sink_interactions['+c+'] :=')
        corig=c.strip('\"')
        for neighbors in big_PPI.predecessors(corig+suff+'_sink'):
            if neighbors in trares:
                file.write('(\"'+neighbors+'\",\"'+sink+'\")')
        file.write(';\n')


    #CAPACITIES
    #capacities are per edge, regardless of commodity
    file.write('\nparam capacity default 1:=\n')
    all_caps={}
    all_direct_caps={}
    ##for hierarchical capacties, keep source and sink the same
    for i in phenres:
        total_cap=0.0
        for c in commodity_source_weights.keys():
            if i in commodity_source_weights[c].keys():
                total_cap+=float(commodity_source_weights[c][i])
        all_caps[i]=total_cap
        #        print 'Source cap to '+i+':'+str(total_cap)
    for i in all_caps.keys():
        file.write(source+' \"'+i+'\"\t'+str(float(all_caps[i]/max(0.000001,sum(all_caps.values()))))+'\n')

    ##now add capacities for direct edges, if they exist...
    for i in directres:
        total_cap=0.0
        for c in commodity_direct_weights.keys():
            if i in commodity_direct_weights[c].keys():
                total_cap+=float(commodity_direct_weights[c][i])
        all_direct_caps[i]=total_cap

    for i in all_direct_caps.keys():
        capsum=sum(all_direct_caps.values())
        file.write(source+' \"'+i+'\"\t'+str(float(all_direct_caps[i]/capsum))+'\n')

    all_tra_caps={}
    for i in trares:
        total_cap=0.0
        for c in commodity_sink_weights.keys():
            if i in commodity_sink_weights[c].keys():
                total_cap+=float(commodity_sink_weights[c][i])
        all_tra_caps[i]=total_cap

#        print 'Sink cap from '+i+':'+str(total_cap)
    for i in all_tra_caps.keys():
        file.write('\"'+i+'\" '+sink+'\t'+str(float(all_tra_caps[i]/sum(all_tra_caps.values())))+'\n')


    ##now add in hierarchical capactities
    other_node_caps={}
    for comm in node_caps.keys():
        for node in node_caps[comm].keys():
            if node in other_node_caps.keys():
                mval=max(other_node_caps[node],node_caps[comm][node])
            else:
                mval=node_caps[comm][node]
            other_node_caps[node]=mval

    ##now write all to file
    for i in other_node_caps.keys():
        if i not in all_tra_caps.keys() and i!=source and i in big_PPI.nodes():
            for suc in big_PPI.successors(i):
                ew=float(big_PPI.get_edge_data(i,suc)['weight'])
                if 'sink' not in suc and suc!=sink and ew!=0.0:
                    file.write('\"'+i+'\" \"'+suc+'\"\t'+str(float(other_node_caps[i]))+'\n')
                
    file.write(';\n')

    #COSTS
    file.write('\nparam cost default 0 :=\n') ##default to 0 to ensure that unselected edges are not used
    #first add costs for source and sink edges, since those vary based on commodity
    #make sure each commodity sums to 1...
    #added weighting here to weight by commodity


    ##normalize al weights to sum of all weights from source
    all_source_weights=float(sum([sum(commodity_source_weights[c].values()) for c in commodity_source_weights.keys()]))
    print 'All source weights '+str(all_source_weights)


    all_sink_weights=float(sum([sum(commodity_sink_weights[c].values()) for c in commodity_sink_weights.keys()]))
    print 'All sink weights '+str(all_sink_weights)

    ##normalize all direct weights separately
    all_direct_weights=float(sum([sum(commodity_direct_weights[c].values()) for c in commodity_direct_weights.keys()]))
    if(all_direct_weights)>0.0:
        print 'All direct source weights '+str(all_direct_weights)

    all_source_nodes=set()
    [all_source_nodes.update(commodity_source_weights[k].keys()) for k in commodity_source_weights.keys()]
#
    for c in commodity_source_weights.keys():
        for node in commodity_source_weights[c]:    
#        for node in all_source_nodes:
#            if node in commodity_source_weights[c].keys():
            file.write('\"'+source+'\" \"'+node+'\" \"'+c+'\" '+str(commodity_source_weights[c][node]/all_source_weights)+'\n')
#            else:
#                print 'Node '+node+' is not connected to Source in commodity '+c+', zeroing'
#                file.write('\"'+source+'\" \"'+node+'\" \"'+c+'\" '+str(0.0)+'\n')
                
    for c in commodity_direct_weights.keys():
        for node in commodity_direct_weights[c]:
            file.write('\"'+source+'\" \"'+node+'\" \"'+c+'\" '+str(commodity_direct_weights[c][node]/all_direct_weights)+'\n')

    for c in commodity_sink_weights.keys():
        for node in commodity_sink_weights[c]:
            if node in trares:
##updated 6/26/14: this penalizes commodities with large amounts of mRNA targets, we don't want to do that.
##                file.write('\"'+node+'\" \"'+sink+'\" \"'+c+'\" '+str(commodity_sink_weights[c][node]/sum(commodity_sink_weights[c].values()))+'\n')
                file.write('\"'+node+'\" \"'+sink+'\" \"'+c+'\" '+str(commodity_sink_weights[c][node]/all_sink_weights)+'\n')

    #the add costs for other edges, simply duplicating values for all commodities
    for two_node_edge in big_PPI.edges():
        n1=two_node_edge[0]
        n2=two_node_edge[1]
        if n1!=source and n2!=sink and n1 not in dummy_nodes and n2 not in dummy_nodes: #make sure we're not a dummy node!
            weight=float(big_PPI.get_edge_data(n1,n2)['weight'])
            if weight>=cap:
                weight=cap
            if weight==0:
#                print 'zero weight1: '+n1+' '+n2
#                weight=0.00001
               continue
            else:
                for c in commodity_names:
                    file.write('\"'+n1+'\" \"'+n2+'\" \"'+c+'\" '+str(weight)+'\n')

    file.write(';')
    file.close()

def write_mcf_amplfile(wholename,solver):
#    print 'Writing ampl file'
    #it is possible to select another solver; we decided to use LOQO but you can use CPLEX etc...
    lines=[]
    lines.append('model '+wholename+'changeflow.mod;\n')
    lines.append('data '+wholename+'.dat;\n')
    lines.append('option solver '+solver+'; \n')
    lines.append('solve; \n')
    lines.append('printf {(i,j) in all_interactions,k in commodities}: "%s\\t%s\\t%s\\t%f\\n", i,j,k,X[i,j,k]>'+wholename+'.txt; \n')
    outFile = open(wholename+'.ampl', 'w')
    for line in lines:
        outFile.write(line)
    outFile.close()
    

#It creates an ampl file that uses the mod and dat file to have the output txt file
def writeamplfile(wholename, solver):
    #it is possible to select another solver; we decided to use LOQO but you can use CPLEX etc...
    lines=[]
    lines.append('model '+wholename+'changeflow.mod;\n')
    lines.append('data '+wholename+'.dat;\n')
    lines.append('option solver '+solver+'; \n')
    lines.append('solve; \n')
    lines.append('printf {(i,j) in interactions}: "%s\\t%s\\t%f\\n", i,j,X[i,j]>'+wholename+'.txt; \n')
    outFile = open(wholename+'.ampl', 'w')
    for line in lines:
        outFile.write(line)
    outFile.close()



# #creates the dat file
# def writeorigdatfile(big_PPI,trares,phenres,mirnares,outputfilename,source,sink, cap,usetargetcapacity=False,diff_ex_vals=dict(),de_cap='sink',node_caps={},debug=False,sinkGamma=False):

#     #print this on the screen, but it is not appearing in the dat file    
#    # print "vertices number:", len(big_PPI.nodes())

#     #open the dat file
#     file=open(outputfilename+'.dat','w')
#     print 'Writing network with '+str(len(big_PPI.nodes()))+' nodes and '+str(len(big_PPI.edges()))+' edges'
#     #print ','.join(big_PPI.nodes())
#     #print source,sink
#     #SUM WEIGHTS SOURCE-PHEN AND MRNA-SINK (USE LATER FOR CAPACITY)
				
#     # sum the weights associated to the genetic and transcriptional data. They will be used for normalization.
#     total_sink_weights=0
#     for i in trares:
#         total_sink_weights=total_sink_weights+big_PPI.get_edge_data(i,sink)['weight']

#     total_source_weights=0
#     for i in phenres:
#         total_source_weights=total_source_weights+big_PPI.get_edge_data(source,i)['weight']


#     #START WRITING .DAT FILE    
#     # start writing on the dat file defining the set of proteins
#     #write all nodes of the graph (proteins and mRNA, but not source and sink)
#     file.write('set proteins := ')
#     proteins=big_PPI.nodes()	
#     for i in proteins:
#         if i!=sink and i!=source:
#             file.write('"'+i+'"'+' ')
#     file.write(';\n')

#     #write down the source and sink
#     file.write('set sourcesink := '+sink+' '+source+' ;\n')



#     #INTERACTIONS
#     #define the interactions
#     file.write('set interactions := ')

#     #since our graph is directed, all we need for interactions is to print information about all edges of the big_PPI (proteins, TF, mRNA, sourcem sink are all included)
		
#     for two_edge_nodes in big_PPI.edges():
#         file.write("(\""+two_edge_nodes[0]+"\",\""+two_edge_nodes[1]+"\")")

#     file.write(';\n')


#     #SOURCE-PHEN AND MRNA-SINK INTERACTIONS
#     file.write('set source_interactions := ')
#     for i in phenres:
#         file.write('('+source+',\"'+i+'\") ')			
#     file.write(';\n')

#     file.write('set sink_interactions :=')
#     for i in trares:
#         file.write('(\"'+i+'\",'+sink+') ')			
#     file.write(';\n')

    
#     #CAPACITIES
#     # write the capacities
#     file.write('\nparam capacity default 1:=\n')
#     for i in phenres:
#         file.write(source+' \"'+i+'\"\t'+str(float(big_PPI.get_edge_data(source,i)['weight'])/(total_source_weights))+'\n')			
#     for i in trares:
#         file.write('\"'+i+'\" '+sink+'\t'+str(float(big_PPI.get_edge_data(i,sink)['weight'])/(total_sink_weights))+'\n')			
#     file.write(';\n')

#     #COST
#     #write the cost
#     #make sure to cap costs
#     file.write('\nparam \tcost :=\n')

    
#     #costs are only between proteins and proteins or proteins and mRNA
#     for protein in big_PPI:
#         if (protein!=source):
#             for neighbors in big_PPI.successors(protein):

#                 if (neighbors!=sink):
#                     weight = float(big_PPI.get_edge_data(protein,neighbors)['weight'])
#                     if weight >= cap:
#                         weight=cap
#                     if weight ==0:
#                         print "zero weight: "+protein +' '+neighbors#
#                         #continue
#                         weight=0.00001

                    
#                     file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(weight)+'\n')
#                 #mrna->sink
#                 else:
#                     file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(float(big_PPI.get_edge_data(protein, neighbors)['weight'])/(total_sink_weights))+'\n')
#         else:
#             #source->prot
#             for neighbors in big_PPI.successors(protein):
#                 file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(float(big_PPI.get_edge_data(protein, neighbors)['weight'])/(total_source_weights))+'\n')

#     file.write(';')
#     file.close()
   

# def writedatfile_with_dummies(big_PPI,trares,phenres,mirnares,outputfilename,source,sink, cap):
#     #print this on the screen, but it is not appearing in the dat file    
#     print "vertices number:", len(big_PPI.nodes())
#     if len(mirnares)==0:
#         'Wrong function, call original writedatfile'
#         return
#     #open the dat file
#     file=open(outputfilename+'.dat','w')

#     #SUM WEIGHTS SOURCE-PHEN AND MRNA-SINK (USE LATER FOR CAPACITY)
				
#     # sum the weights associated to the genetic and transcriptional data. They will be used for normalization.
#     total_sink_weights=0
#     for i in trares:
#         total_sink_weights=total_sink_weights+big_PPI.get_edge_data(i,sink)['weight']

#     total_source_weights=1
    
#     total_indirect_weights=0
#     for i in phenres:
#         total_indirect_weights=total_indirect_weights+big_PPI.get_edge_data('indirect',i)['weight']
    
#     total_direct_weights=0
#     if len(mirnares)>0:
#         for i in mirnares:
#             total_direct_weights=total_direct_weights+big_PPI.get_edge_data('direct',i)['weight']
            
#     print 'sink weights:',str(total_sink_weights),'indirect weights:',str(total_indirect_weights),'direct weights:',str(total_direct_weights)

#     #START WRITING .DAT FILE    
#     # start writing on the dat file defining the set of proteins
#     #write all nodes of the graph (proteins and mRNA, but NOT source and sink)
#     file.write('set proteins := ')
#     proteins=big_PPI.nodes()	
#     for i in proteins:
#         if i!=sink and i!=source:
#             file.write('\"'+i+'\" ')
#     file.write(';\n')

#     #write down the source and sink
#     file.write('set sourcesink := '+sink+' '+source+' ;\n')


#     #INTERACTIONS
#     #define the interactions
#     file.write('set interactions := ')

#     #since our graph is directed, all we need for interactions is to print information about all edges of the big_PPI (proteins, TF, mRNA, sourcem sink are all included)
		
#     for two_edge_nodes in big_PPI.edges():
#         file.write("(\""+two_edge_nodes[0]+"\",\""+two_edge_nodes[1]+"\")")

#     file.write(';\n')


#     #SOURCE-PHEN AND MRNA-SINK INTERACTIONS
#     file.write('set source_interactions := ')

#     ##if we're using mirnas, connect source to direct/indirect. 
#     file.write('('+source+','+'direct) ')
#     file.write('('+source+','+'indirect) ')

#     file.write(';\n')

#     file.write('set sink_interactions :=')
#     for i in trares:
#         file.write('(\"'+i+'\",'+sink+') ')			
#     file.write(';\n')

    
#     #CAPACITIES
#     # write the capacities -- only for edges going out of source or into sink
#     file.write('\nparam capacity default 1:=\n')
#     file.write(source+' '+'direct'+'\t'+'0.5'+'\n')
#     file.write(source+' '+'indirect'+'\t'+'0.5'+'\n')
    
#     for i in trares:
#         file.write('\"'+i+'\" '+sink+'\t'+str(float(big_PPI.get_edge_data(i,sink)['weight'])/(total_sink_weights))+'\n')


    
#     ##now add in hierarchical capactities
#     other_node_caps={}
#     for comm in node_caps.keys():
#         for node in node_caps[comm].keys():
#             if node in other_node_caps.keys():
#                 mval=max(other_node_caps[node],node_caps[comm][node])
#             else:
#                 mval=node_caps[comm][node]

#     ##now write all to file
#     for i in other_node_caps.keys():
#         if i not in all_tra_caps.keys() and i!=source:
#             for suc in networkx.neighbors(i):
#                 file.write('\"'+i+'\" '+suc+'\t'+str(float(other_node_caps[i]))+'\n')
                
#     file.write(';\n')

#     #COST
#     #write the cost
#     #make sure to cap costs
#     file.write('\nparam:\tcost :=\n')

# ##    if(len(mirnares)>0):
#  ##       for i in phenres:
#   ##          file.write('indirect'+' '+i+'\t'+str(float(big_PPI.get_edge_data('indirect',i)['weight'])/(total_indirect_weights))+'\n')	
# #        for i in mirnares:
# #            file.write('direct'+' '+i+'\t'+str(float(big_PPI.get_edge_data('direct',i)['weight'])/(total_direct_weights))+'\n')

#     #costs are only between proteins and proteins or proteins and mRNA
#     for protein in big_PPI:
#         if (protein!=source):
#             for neighbors in big_PPI.successors(protein):
#                 if (neighbors!=sink):
#                     weight = float(big_PPI.get_edge_data(protein,neighbors)['weight'])
#                     if weight >= cap:
#                         weight=cap
#                     if weight ==0:
#                         print "zero weight: "+protein +' '+neighbors
#                         weight=0.00001
#                     file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(weight)+'\n')
#                 #mrna->sink
#                 else:
#                     file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(float(big_PPI.get_edge_data(protein, neighbors)['weight'])/(total_sink_weights))+'\n')
#         elif protein==source:
#             #source->prot
#             for neighbors in big_PPI.successors(protein):
#                 file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(float(big_PPI.get_edge_data(protein, neighbors)['weight'])/(total_source_weights))+'\n')
# #        elif protein=='direct' and len(mirnares)>0:
# #            for neighbors in big_PPI.neighbors(protein):
# #                file.write(protein+' '+neighbors+'\t'+str(float(big_PPI.get_edge_data(protein,neighbors)['weight'])/(total_indirect_weights))+'\n')
# #        else:
# #            for neighbors in big_PPI.neighbors(protein):
# #                file.write(source+' '+neighbors+'\t'+str(float(big_PPI.get_edge_data(protein,neighbors)['weight'])/(total_indirect_weights))+'\n')

#     file.write(';')
#     file.close()


# def writedatfile_with_multiple_treatments(big_PPI,trares,phenres,mirnares,outputfilename,source,sink,cap,usetargetcapacity=True,diff_ex_vals=dict(),de_cap='sink',node_caps={},debug=False):
#     '''
#     Third version of writedatafile_* that incorporates multiple treatments
#     '''


#     #print this on the screen, but it is not appearing in the dat file    
#     print "vertices number:", len(big_PPI.nodes())
#     #open the dat file

#     file=open(outputfilename+'.dat','w')

#     #print 'Source:'+source,', sink:'+sink
#     #SUM WEIGHTS SOURCE-PHEN AND MRNA-SINK (USE LATER FOR CAPACITY)
				
#     # sum the weights associated to the genetic and transcriptional data. They will be used for normalization.
#     total_sink_weights=0
#     for treat in big_PPI.predecessors(sink):
#         total_sink_weights=total_sink_weights+big_PPI.get_edge_data(treat,sink)['weight']

#     total_source_weights=0 ##
#     for treat in big_PPI.successors(source):
#         total_source_weights=total_source_weights+big_PPI.get_edge_data(source,treat)['weight']

            
#     print 'sink weights:',str(total_sink_weights),'source weights: '+str(total_source_weights)#indirect weights:',str(total_indirect_weights),'direct weights:',str(total_direct_weights)

#     #START WRITING .DAT FILE    
#     # start writing on the dat file defining the set of proteins
#     #write all nodes of the graph (proteins and mRNA, but NOT source and sink)
#     file.write('set proteins := ')
#     proteins=big_PPI.nodes()	
#     for i in proteins:
#         if i!=sink and i!=source:
#             file.write('\"'+i+'\" ')
#     file.write(';\n')

#     #write down the source and sink
# #    file.write('set sourcesink := '+sink+' '+source+' ;\n')
#     file.write('set sourcesink := \"'+sink+'\" \"'+source+'\" ')
#     if usetargetcapacity:
#         all_neighbors=big_PPI.successors(source)
#         all_neighbors.update(big_PPI.predecessors(sink))
#         for treat in all_neighbors:
#             file.write('\"'+treat+'\" ')
# #        for treat_s in big_PPI.predecessors(sink):
# #            file.write('\"'+treat_s+'\" ')
        
#     file.write(';\n')


#     #INTERACTIONS
#     #define the interactions
#     file.write('set interactions := ')

#     #since our graph is directed, all we need for interactions is to print information about all edges of the big_PPI (proteins, TF, mRNA, sourcem sink are all included)
		
#     for two_edge_nodes in big_PPI.edges():
#         file.write("(\""+two_edge_nodes[0]+"\",\""+two_edge_nodes[1]+"\")")

#     file.write(';\n')


#     #SOURCE-PHEN AND MRNA-SINK INTERACTIONS
#     file.write('set source_interactions := ')

#     ##if we're using mirnas, connect source to direct/indirect. 
#     for treat in big_PPI.successors(source):
#         file.write('(\"'+source+'\",\"'+treat+'\") ')
#         if usetargetcapacity:
#             for targ in big_PPI.neighbors(treat):
#                 file.write('(\"'+treat+'\",\"'+targ+'\") ')


#     file.write(';\n')

#     file.write('set sink_interactions :=')
#     for treat in big_PPI.predecessors(sink):
#         file.write('(\"'+treat+'\",\"'+sink+'\") ')
#         if usetargetcapacity:
#             for i in big_PPI.predecessors(treat):
#                 file.write('(\"'+i+'\",'+treat+') ')			
#     file.write(';\n')

    
#     #CAPACITIES
#     # write the capacities -- only for edges going out of source or into sink
#     file.write('\nparam capacity default 1:=\n')

#     missed_caps=0
#     total_caps=0
#     other_node_caps={}
#     for treat in big_PPI.successors(source):#phenres.keys():
#         #if usetargetcapacity:
#         #    for targ in big_PPI.neighbors(treat):
#         #        file.write('\"'+treat+'\" \"'+targ+'\"\t'+str(float(big_PPI.get_edge_data(treat,targ)['weight']))+'\n')
#         ##now add in hierarchical capactities
#         ##now write all to file

                    
#         ##this is oana's code, not using right now###############################
#         #don't need this anymore...
#         if de_cap=='none' or de_cap=='sink':
#             file.write('\"'+source+'\" \"'+treat+'\"\t'+str(float(big_PPI.get_edge_data(source,treat)['weight'])/(total_source_weights))+'\n')
#         elif de_cap=='source':
#             #first check to see if the diff_ex_values is in the keys, this means that we need to select
#             if source in diff_ex_vals.keys():##single source
#                 dev=diff_ex_vals[source]
#                 #print 'Got '+str(len(dev))+' diff ex vals for '+source
#                 if treat in dev.keys():
#                     val=dev[treat]/sum(dev.values())
#                 elif treat not in phenres+mirnares:
#                     print treat+' not in network'
#                 else:
#                     val=float(sum(dev.values()))/float(len(dev.values()))
#                     missed_caps+=1
#                     if(debug):
#                         print 'No differential expression value for '+treat+', using '+str(val)     
#                 total_caps+=1
#                 file.write('\"'+source+'\" \"'+treat+'"\t'+str(val)+'\n')
#             elif treat in diff_ex_vals.keys():##multiple sources
#                 dev=diff_ex_vals[treat]
#                 print 'Got '+str(len(dev))+' diff ex vals for '+treat
#                 for targ in big_PPI.successors(treat):
#                     if targ in dev.keys():
#                         val=dev[targ]/sum(dev.values())
#                     elif targ not in phenres+mirnares:
#                         print targ+' not in network'
#                     else:
#                         val=float(sum(dev.values()))/float(len(dev.values()))
#                         missed_caps+=1
#                         if debug:
#                             print 'No differential expression value for '+targ+', using '+str(val)
#                     total_caps+=1
#                     file.write('\"'+treat+'\" \"'+targ+'"\t'+str(val)+'\n')
#             else:
#                 "No differential expression values for node "+source+" or "+treat


#         else:
#             print ("Diffex capacities for all nodes not yet implemented")
        
#     #End diffex capping stuff#############################################
#     if debug:
#         print 'Missed '+str(missed_caps)+' diffex vals for '+source+' out of '+str(total_caps)

#     ##now write out node capacities, only for 1 source
#     for node in node_caps[source].keys():
# #        if node in other_node_caps.keys():
# #            mval=max(other_node_caps[node],node_caps[source][node])
# #        else:
#         other_node_caps[node]=node_caps[source][node]
#         #print 'Node capacity for '+node+' in commodity '+source+' is '+str(mval)

#     for i in other_node_caps.keys():
#         if i not in big_PPI.predecessors(sink) and i!=source:
#             for suc in big_PPI.successors(i):
#                 #print 'Edge capacity from '+i+' to '+suc+' is '+str(float(other_node_caps[i]))
#                 file.write('\"'+i+'\" \"'+suc+'\"\t'+str(float(other_node_caps[i]))+'\n')

#     for treat in big_PPI.predecessors(sink):
#         file.write('\"'+treat+'\" \"'+sink+'\"\t'+str(float(big_PPI.get_edge_data(treat,sink)['weight'])/(total_sink_weights))+'\n')	
#         if usetargetcapacity:
#             for targ in big_PPI.predecessors(treat):
#                 file.write('\"'+targ+'\" \"'+treat+'\"\t'+str(float(big_PPI.get_edge_data(targ,treat)['weight']))+'\n')
#     file.write(';\n')

#     #COST
#     #write the cost
#     #make sure to cap costs
#     file.write('\nparam:\tcost :=\n')

#     #costs are only between proteins and proteins or proteins and mRNA
#     for protein in big_PPI:
#         if (protein!=source):
#             for neighbors in big_PPI.successors(protein):
#                 if (neighbors!=sink):
#                     weight = float(big_PPI.get_edge_data(protein,neighbors)['weight'])
#                     if weight >= cap:
#                         print 'Capped ('+str(cap)+') weight: '+protein+' '+neighbors+' '+str(weight)
#                         weight=cap

#                     if weight ==0:
#                         weight=0.0001
#                         print "zero weight: "+protein +' '+neighbors
#                     file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(weight)+'\n')
#                 #mrna->sink
#                 else:
#                     file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(float(big_PPI.get_edge_data(protein, neighbors)['weight'])/(total_sink_weights))+'\n')
#         elif protein==source:
#             #source->prot
#             for neighbors in big_PPI.successors(protein):
#                 file.write('\"'+protein+'\" \"'+neighbors+'\"\t'+str(float(big_PPI.get_edge_data(protein, neighbors)['weight'])/(total_source_weights))+'\n')

#     file.write(';')
#     file.close()


#def write_all_files(PPI_with_weights,trares,phenres,directres,output,source,sink,cap,gamma,solver,usetargetcapa# city=False,diff_ex_vals=dict(),de_cap='sink',node_caps={},debug=False):
# def write_all_files(PPI_with_weights,trares,phenres,output,source,sink,cap,gamma,solver,usetargetcapacity=False,diff_ex_vals=dict(),de_cap='sink',node_caps={},debug=False,sinkGamma=False):
# #    print de_cap
#     writeorigdatfile(PPI_with_weights, trares, phenres,output,source,sink,cap,usetargetcapacity,diff_ex_vals,de_cap,node_caps,debug)        
# #    writedatfile_with_multiple_treatments(PPI_with_weights, trares, phenres,directres,output,source,sink,cap,usetargetcapacity,diff_ex_vals,de_cap,node_caps,debug)        

#     writechangeflow(output,gamma,sinkGamma)
#     # create ampl file
#     writeamplfile(output,solver)

#def write_mcf_files(PPI_with_weights,trares,phenres,directres,output,source,sink,cap,gamma,solver,usetargetcapacity=False,diff_ex_vals=dict(),de_cap='sink',node_caps={},debug=False):
def write_mcf_files(PPI_with_weights,trares,phenres,output,source,sink,cap,gamma,solver,usetargetcapacity=False,diff_ex_vals=dict(),de_cap='sink',node_caps={},debug=False,sinkGamma=False):
#    write_mcf_datfile(PPI_with_weights,trares,phenres,directres,output,source,sink,cap,usetargetcapacity,diff_ex_vals,de_cap,node_caps,debug)
    write_mcf_datfile(PPI_with_weights,trares,phenres,output,source,sink,cap,usetargetcapacity,diff_ex_vals,de_cap,node_caps,debug)
    write_mcf_changeflow(output,gamma,sinkGamma)
    write_mcf_amplfile(output,solver)

