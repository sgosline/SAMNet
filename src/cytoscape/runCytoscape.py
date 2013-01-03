'''
runCytoscape.py
Helper file to facilitate plotting of SAMNet output in cytoscape


Copyright (c) 2012 Sara JC Gosline
sgosline@mit.edu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

'''


import os
import sys
from optparse import OptionParser

def main():
    '''
    handles various types of identifier matching
    '''
    parser=OptionParser()
    parser.add_option('--nomrna',default=False,action='store_true',dest='nomrna',help='Set if you want no mRNA in the network')
    parser.add_option('--has_DE_info',default=False,action='store_true',dest='has_de',help='Set this argument if you want to use differential expression info')
    parser.add_option('--isMCF',default=False,action='store_true',dest='mcf',help='Set if you want multi-commodity view of network')
    parser.add_option('--noSymbol',default=False,action='store_true',dest='nosymbol',help='Set this flag to use original sifs, not symbol-replaced ones')

    (opts,args)=parser.parse_args()

    if len(args)<1:
        exit("Need argument for network file prefix")

    
    prefix=args[0]

    fpath=os.path.dirname(os.path.abspath( __file__ ))

    if opts.nosymbol:
        sym=''
    else:
        sym='_symbol'
        
    eattr=os.path.realpath(prefix+'_ppi_attributes'+sym+'.eda')        
    propFile=fpath+'/nodeEdgeVizMapFile.props'
    if opts.nomrna:
        net=os.path.realpath(prefix+'_no_mrna'+sym+'.sif')

    elif opts.mcf:
        net=os.path.realpath(prefix+'_mcfs'+sym+'.sif')
        eattr=os.path.realpath(prefix+'_edge_commodity'+sym+'.eda')
        ##extra fiels for mcf
        eattr2=os.path.realpath(prefix+'_edge_type'+sym+'.eda')
        node_tab=os.path.realpath(prefix+'_node_comm_flow'+sym+'.noa')
        propFile=fpath+'/nodeEdgeVizMapFile_mcf.props'
    else:
        net=os.path.realpath(prefix+'_all'+sym+'.sif')


    nattr1=os.path.realpath(prefix+'_node_type'+sym+'.noa')
    nattr2=os.path.realpath(prefix+'_node_flow'+sym+'.noa')



    if opts.has_de:
        nattr3=os.path.realpath(prefix+'_DiffExpr.noa')
        cmd='Cytoscape -N '+net+' -n '+nattr1+' -n '+nattr2+' -n '+nattr3+' -e '+eattr+' -V '+propFile
    ##now check to see if sif file exists
    else:
        cmd='Cytoscape -N '+net+' -n '+nattr1+' -n '+nattr2+' -e '+eattr+' -V '+propFile

    if opts.mcf:
        cmd+=' -m '+node_tab+' -e '+eattr2

    if(os.path.exists(net)):
        os.system(cmd)
    else:
        print "File",net,"does not exist"
       
   

if __name__ == "__main__" :
    main()
