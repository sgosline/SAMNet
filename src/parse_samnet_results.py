#!/usr/local/python2.6/bin/python

'''
Parse the result files created by the SAMNet into a viewable web page

'''
import os
import sys
#import networkx
#import Shared.LoggerFactory
#import Utils.XGraph2SIF
#import CFA3Config
#import CFA3Util
#import Constants
#import PCSTResults
#import make_pcst_input
#import pickle
from collections import defaultdict

__author__ = 'Shao-shan Carol Huang,Sara Gosline'
__email__ = 'shhuang@mit.edu,sgosline@mit.edu'


EBE_CLUST_ITER = 10
def title_html(outdir):
    '''
    This creates the title frame, pretty basic so far
    '''
    file = open(outdir+'/title.html','w')
    file.writelines("""<center><h2><a href="http://fraenkel.mit.edu/samnetweb/" target="_parent">S A M N E T</a></h1></center>""")
    file.writelines("""<center><h3>A multi-commodity flow-based data integration tool</h3></center>""")
    file.close()

def summary_html(outdir):
    '''
    This creates the summary panel to the right with all the run-specific details
    '''
    ##unique sets of protein and mrna weights
    terminal = set([arr.strip().split('\t')[1] for arr in open(outdir+'/../data/proteinWeights','rU').readlines()])
    mrnaterminal = set([arr.strip().split('\t')[1] for arr in open(outdir+'/../data/exp','rU').readlines()])

    ##get original NONSYMBOL file
    sif_file = outdir + "/samnetoutmultiComm_mcfs.sif"
    file = open(sif_file)
    tindex = 0
    mrnatindex = 0
    nontindex = 0
    nodeSet = set()
    while 1:
        line = file.readline()
        if line == "": break
        temp = line.strip().split()
        nodeSet.add(temp[0])
        nodeSet.add(temp[2])
    file.close()
    for node in nodeSet:
        if node in terminal:
            tindex += 1
        if node in mrnaterminal:
            mrnatindex += 1
    ##count how many of original results are included
    mrnanontindex = len(mrnaterminal)-mrnatindex
    nontindex = len(terminal) - tindex

    #this page is the cytoscape panel
    file = open(outdir+'/CytoscapeWebControls.html','w')
    file.writelines("""
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">\n
    <html><head>
    <title>Cytoscape Controls</title>



    <script type="text/javascript" src="http://www.google.com/jsapi"></script>
    <script type="text/javascript">
      google.load('visualization', '1', {packages: ['corechart']});
    </script>
    <script type="text/javascript">
      function drawVisualization() {
        // Create and populate the data table.
        var data = new google.visualization.DataTable();
        data.addColumn('string', 'Task');
        data.addColumn('number', 'Included Source Nodes');
        data.addRows(5);
        data.setValue(0, 0, 'Included Source (e.g. Protein Weight) Nodes');
        data.setValue(0, 1, %d);
        data.setValue(1, 0, 'Excluded Source (e.g. Protein Weight) Nodes');
        data.setValue(1, 1, %d);
      
        // Create and draw the visualization.
        new google.visualization.PieChart(document.getElementById('visualization')).
            draw(data, {title:"Protein Weight Nodes"});
      }
      google.setOnLoadCallback(drawVisualization);
    </script>

    <script type="text/javascript">
      function drawVisualizationDna() {
        // Create and populate the data table.
        var data = new google.visualization.DataTable();
        data.addColumn('string', 'Task');
        data.addColumn('number', 'Included Terminals');
        data.addRows(5);
        data.setValue(0, 0, 'Included Sink (e.g. mRNA) Nodes');
        data.setValue(0, 1, %d);
        data.setValue(1, 0, 'Excluded Sink (e.g. mRNA) Nodes');
        data.setValue(1, 1, %d);

        // Create and draw the visualization.
        new google.visualization.PieChart(document.getElementById('visualizationdna')).
            draw(data, {title:"Sink Nodes"});
      }
      google.setOnLoadCallback(drawVisualizationDna);
    </script>


    <script type="text/javascript"> 
    function fullscreen() { 
    winRef = window.open("samnetmultiComm_result.html", "CytoscapeFullscreenWindow");
    }
    </script>

    <script type="text/javascript">
    function selectForceDirectedLayout() {
        document.getElementById('forceDirectedLayout').style.fontWeight='bold';
        document.getElementById('radialLayout').style.fontWeight='normal';
        document.getElementById('circleLayout').style.fontWeight='normal';
        vis.layout('ForceDirected');
    }
    </script>
</head>
    <body >
<center><h3>Legend for the Network Visualization</h3><img src="../../../Legend.png" height="150" width="225"/><br></center>

<table>
<tr><td><h3>Input Data</h3>
<tr><td>
    <a href="../data/proteinWeights" >Protein weights inputs</a><br/>
    <a href="../data/exp" >mRNA expression inputs</a><br/>
    <a href="../data/tf2gene" >Transcription Factor to DNA interaction file</a>
    <br/><a href="../data/interactome" >Protein-Protein Interaction PKL file</a><br/>
</td></tr>
<tr><td><h3>Output Files</h3>
    <a href="samnetoutmultiComm_mcfs_symbol.sif">Download SAMNet network in SIF format (Cytoscape) </a><br/>
    <a href="samnetoutmultiComm_edge_commodity_symbol.eda">Download Edge Flow values </a><br/>
    <a href="samnetoutmultiComm_edge_type_symbol.eda">Download Edge Type values </a><br/>
    <a href="samnetoutmultiComm_node_flow_symbol.noa">Download Node Flow values </a><br/>
    <a href="samnetoutmultiComm_node_type_symbol.noa">Download Node Type values </a><br/>
    <a href="../samnet_CYTOEDA_Benjamini_1.0_sigDavidTerms.xls">Download DAVID terms for each commodity</a><br/></td></tr>
</table>
    
    <h3>Input nodes in the Solution</h3><center>
<table><tr><td>
  <div id="visualization" style="width: 300px; height: 200px;"></div>
</td>
</tr>
<tr><td>
  <div id="visualizationdna" style="width: 300px; height: 200px;"></div>
</td></tr>
</table>
    </body></html>
""" %(tindex,nontindex,mrnatindex,mrnanontindex))

    file.close()

##    <table><tr><td><img src="http://chart.apis.google.com/chart?chs=400x200&cht=p3&chd=t:%d,%d&chdl=Terminals+Included|Terminals+Excluded&chdlp=b&chp=0.4&chtt=Terminal+Nodes" width="200" height="100" alt="Terminal Nodes" /></td><td><img src="http://chart.apis.google.com/chart?chs=400x200&cht=p3&chd=t:%d,%d&chdl=DNA+Terminals+Included|DNA+Terminals+Excluded&chdlp=b&chp=0.4&chtt=DNA+Terminal+Nodes" width="200" height="100" alt="DNA Terminal Nodes" /></td></tr></center></table>

'''
SAMNet post-processing does a lot less than PCST, let's try to remove those that are not necessary
'''
def get_expr_filename(outdir):
    return outdir + "/../data/exp"
def get_protweight_filename(outdir):
    return outdir + "/../data/proteinWeights"
def get_tf_filename(outdir):
    return outdir + "/../data/tf2gene"

def get_html_result(outdir,input_filename):
    return os.path.join(outdir, input_filename+'_result.html')
def get_html_filename(outdir,input_filename):
    return os.path.join(outdir,input_filename+'.html')


##here are the basic files we need
symbol='_symbol' ##set to '' if there are no symbols

#input_filename='samnetoutmultiComm'
def get_sif_filename(outdir,input_filename):
    return os.path.join(outdir,input_filename+'_mcfs'+symbol+'.sif')

def get_node_flow_noa_filename(outdir,input_filename):
    return os.path.join(outdir,input_filename+'_node_flow'+symbol+'.noa')

def get_node_type_noa_filename(outdir,input_filename):
    return os.path.join(outdir,input_filename+'_node_type'+symbol+'.noa')

def get_edge_commodity_eda_filename(outdir,input_filename):
    return os.path.join(outdir,input_filename+'_edge_commodity'+symbol+'.eda')

def get_edge_type_eda_filename(outdir,input_filename):
    return os.path.join(outdir,input_filename+'_edge_type'+symbol+'.eda')


def networkPrep(nodeDict, edgeDict):
    '''
    This function takes the node type dictionary and the edge-commodity dictionary and
    moves it into a cytoscape network string
    '''
    stringstart = "                var networ_json = {\n                    "
    stringstart+="""dataSchema : {\n                        nodes: [{name : "type", type: "string"},\n                               {name : "flow", type: "string"}],\n                   edges: [{name: "interaction",type: "string"}]},\n                           """
    stringstart+="""data: {\n               nodes: [ """

    for node in nodeDict.keys():
    
        ntype = nodeDict[node]['type']
        nflow = nodeDict[node]['flow'] 
        
        node = """{ id: "%s", type: "%s", flow: "%s"},""" % (node, ntype, nflow)
        stringstart += node

    stringstart = stringstart[:-1] + "],\n                         edges: ["
    for comm in edgeDict.keys():
        for edge in edgeDict[comm]:
            e1 = edge[0]
            e2 = edge[1]
            
            e = """{ target: "%s", source: "%s", interaction: "%s" },""" % (e1,e2,comm) #(edge[0], edge[1])
            stringstart += e
    stringstart = stringstart[:-1] + """]\n}\n};"""

    return stringstart

def generateColors(commodities):
    dozencolors=["#003366","#339933","#FFFF66","#FF9999","#9999FF","#66FFFF","#FF9933","#990000","#595959","#00FFFF","#CC9900","#996633"]
    commdict={}
    for ind,c in enumerate(commodities):
        commdict[c]=dozencolors[ind % len(dozencolors)]
    return commdict
    

def sifParser(outdir, input_filename):
    nodeset = set()
    edgedict=defaultdict(list)
    ##i dont think i need these three files
    pwfile = get_protweight_filename(outdir)
    exprfile = get_expr_filename(outdir)
    tf_file = get_tf_filename(outdir)
    
    ##this gets the symbol filename
    siffile = get_sif_filename(outdir,input_filename)
    
    ##but we really need the edges and nodes by commodity
    #edge_comm_filename=get_edge_commodity_eda_filename(outdir,input_filename)
    node_type_filename=get_node_type_noa_filename(outdir,input_filename)
    node_flow_filename=get_node_flow_noa_filename(outdir,input_filename)

    ##first get edges, interaction types
    file = open(siffile).readlines()
    for row in file:
        p1,comm,p2=row.strip().split()

        if p1 in ['S1','T1'] or p2 in  ['S1','T1']:##lets not show the source or sink
            continue
        edgedict[comm].append([p1,p2])
        nodeset.add(p1)
        nodeset.add(p2)
    ##now get node type/flow dictionaries
    nodet,nodef={},{}
    for row in open(node_type_filename,'rU').readlines():
        arr=row.strip().split()
        if len(arr)!=3:
            continue
        nodet[arr[0]]=arr[2]
    for row in open(node_flow_filename,'rU').readlines():
        arr=row.strip().split()
        if len(arr)!=3:
            continue
        nodef[arr[0]]=arr[2]
    nodedict={}
    for node in nodeset:
        t=''
        f=''
        if node in nodet.keys():
            t=nodet[node]
        if node in nodef.keys():
            f=nodef[node]
        nodedict[node]={'type':t,'flow':f}

    return networkPrep(nodedict,edgedict), edgedict.keys()


def result_html_prepare(outdir):
    '''
    This is the main file, with the 3 inset frames
    '''
    input_filename=os.path.join(outdir,'samnetoutmultiComm')

    resultfilename = get_html_result(outdir, input_filename)

    htmlfilename = os.path.basename(get_html_filename(outdir,input_filename))

    resultfile = open(resultfilename,'w')
    resultfile.writelines("""
    <html>
    <title>SAMNetWeb - Results</title>
    <frameset rows="10%,90%">
    <frame src="title.html" name="north" scrolling="no"/>
    <frameset cols="70%,30%">
    """)
    ##now add in html filename
    resultfile.write('<frame src="'+htmlfilename+'" name="center"/>\n')
    ##now add in the rest
    resultfile.writelines("""
    <frame src="CytoscapeWebControls.html" name="east"/>
    </frameset>
    
    </frameset>
    </html>

    """)
    resultfile.close()


def html_prepare(outdir):
    '''
    This is the cytoscape plugin page
    '''
    input_filename=os.path.join(outdir,'samnetoutmultiComm')

    htmlfilename = get_html_filename(outdir,input_filename)


    htmlfile = open(htmlfilename,'w')
    htmlfile.writelines("""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">\n""")
    htmlfile.writelines("""<html>\n""")

    htmlfile.writelines("""    <head>\n""")
    htmlfile.writelines("""        <title>SAMNet Web Output Page</title>\n""")
#    htmlfile.writelines("""        <div id="cytoscapeweb" style="width:100%; height:100%;">\n""")
#    htmlfile.writelines("""            Cytoscape Web will replace the contents of this div with your graph.\n""")
#    htmlfile.writelines("""        </div>\n""")

    htmlfile.writelines("""        <!-- JSON support for IE (needed to use JS API) --> \n""")
    htmlfile.writelines("""        <script type="text/javascript" src="../../../cyto/js/min/json2.min.js"></script>\n""")

    htmlfile.writelines("""        <!-- Flash embedding utility (needed to embed Cytoscape Web) -->\n""")
    htmlfile.writelines("""        <script type="text/javascript" src="../../../cyto/js/min/AC_OETags.min.js"></script>\n""")

    htmlfile.writelines("""        <!-- Cytoscape Web JS API (needed to reference org.cytoscapeweb.Visualization) -->\n""")
    htmlfile.writelines("""        <script type="text/javascript" src="../../../cyto/js/min/cytoscapeweb.min.js"></script>\n""")

    htmlfile.writelines("""        <script type="text/javascript">\n""")
    htmlfile.writelines("""            window.onload=function() {\n""")
    htmlfile.writelines("""                // id of Cytoscape Web container div\n""")
    htmlfile.writelines("""                var div_id = "cytoscapeweb";\n""")

    htmlfile.writelines("""                // you could also use other formats (e.g. GraphML) or grab the network data via AJAX\n""")

    networkline,commodities = sifParser(outdir, input_filename)
    htmlfile.writelines(networkline)
    comcolors=generateColors(commodities)
    ##this describes the node types: source, protein, tf, mRNA
    htmlfile.writelines("""                var shapeMapper = {\n""")
    htmlfile.writelines("""                    attrName: "type",\n""")
    htmlfile.writelines("""                    entries: [ {attrValue: "mrna", value: "DIAMOND"},\n""")
    htmlfile.writelines("""                               {attrValue: "phenotypic", value: "RECTANGLE"},\n""")
    htmlfile.writelines("""                               {attrValue: "mrna,phenotypic", value: "RECTANGLE"},\n""")
    htmlfile.writelines("""                               {attrValue: "transcriptionfactorphenotypic", value: "RECTANGLE"},\n""")
    htmlfile.writelines("""                               {attrValue: "transcriptionfactor", value: "TRIANGLE"},\n""")
    htmlfile.writelines("""                               {attrValue: "transcriptionfactor,mrna", value: "TRIANGLE"}\n""")
    htmlfile.writelines(""" ]\n""")
    htmlfile.writelines("""                };\n""")
    
    ##this describes the node colors
    ##TODO add edge color by commodity
    htmlfile.writelines("""                var colorMapper = {\n""")
    htmlfile.writelines("""                    attrName: "interaction",\n""")
    htmlfile.writelines("""                    entries: [ {attrValue: "", value: "#00CED1"}""")
    for comm in comcolors.keys():
        htmlfile.write(',\n                               {attrValue: "'+comm+'", value: "'+comcolors[comm]+'"}')
#    htmlfile.writelines("""                               {attrValue: "tgfb", value: "#E8E8E8"},\n""")
#    htmlfile.writelines("""                               {attrValue: "tgfb", value: "#E8E8E8"},\n""")
    htmlfile.writelines(""" ]\n""")
    htmlfile.writelines("""                };\n""")
    htmlfile.writelines("""                var visual_style = {\n""")
    htmlfile.writelines("""                    nodes: {\n""")
#    htmlfile.writelines("""                        color: { discreteMapper: colorMapper },\n""")
    htmlfile.writelines("""                        shape: { discreteMapper: shapeMapper },\n""")
    htmlfile.writelines("""                        label: { passthroughMapper: {attrName: "id"}},\n""")
    htmlfile.writelines("""                        borderWidth: 3,\n""")
    htmlfile.writelines("""                        size: 30,\n""")
    htmlfile.writelines("""                        labelVerticalAnchor: "bottom",\n""")
    htmlfile.writelines("""                        labelFontSize: 16\n""")

    htmlfile.writelines("""                     },\n""")
    htmlfile.writelines("""                     edges: {\n""")
    htmlfile.writelines("""                         color: { discreteMapper: colorMapper },\n""")
    htmlfile.writelines("""                         width: 2\n""")
    htmlfile.writelines("""                     }\n""")

    htmlfile.writelines("""                };\n""")

    htmlfile.writelines("""                var layout = {\n""")
    htmlfile.writelines("""                    name: "ForceDirected"\n""")

    htmlfile.writelines("""                };\n""")

    htmlfile.writelines("""                // initialization options\n""")

    htmlfile.writelines("""                var options = {\n""")
    htmlfile.writelines("""                    // where you have the Cytoscape Web SWF\n""")
    htmlfile.writelines("""                    swfPath: "../../../cyto/swf/CytoscapeWeb",\n""")
    htmlfile.writelines("""                    // where you have the Flash installer SWF\n""")
    htmlfile.writelines("""                    flashInstallerPath: "../../../cyto/swf/playerProductInstall"\n""")
    htmlfile.writelines("""                };\n""")

    htmlfile.writelines("""                // init and draw\n""")
    htmlfile.writelines("""                var vis = new org.cytoscapeweb.Visualization(div_id, options);\n""")
    htmlfile.writelines("""//                vis.ready(function() {\n""")
    htmlfile.writelines("""//                    document.getElementById("layout").onclick = function() {\n""")
    htmlfile.writelines("""//                        vis.layout("CompoundSpringEmbedder");\n""")
    htmlfile.writelines("""//                    };\n""")
    htmlfile.writelines("""//                });\n""")
    htmlfile.writelines("""                var draw_options = {\n""")
    htmlfile.writelines("""                    network: networ_json,\n""")

    htmlfile.writelines("""                    layout: "ForceDirected",\n""")
    htmlfile.writelines("""                    visualStyle: visual_style,\n""")
    htmlfile.writelines("""                    //panZoomControlVisible: false\n""")
    htmlfile.writelines("""                };\n""")
    htmlfile.writelines("""                vis.draw(draw_options);\n""")
    htmlfile.writelines("""            };\n""")
    htmlfile.writelines("""            graphResize()\n""")
    htmlfile.writelines("""            showhide('cytoscapeweb');\n""")

    htmlfile.writelines("""        </script>\n""")

    htmlfile.writelines("""    </head>\n""")

    htmlfile.writelines("""    <body>\n""")
#    htmlfile.writelines("""        <div id="cytoscapeweb" style="width:100%; height:100%;">\n""")
    htmlfile.writelines("""        <div id="cytoscapeweb" style="position:absolute; left:0px; top:0px; width:100%; height:100%; z-index:1;">\n""")
    htmlfile.writelines("""            Cytoscape Web will replace the contents of this div with your graph.\n""")
    htmlfile.writelines("""        </div>\n""")
    htmlfile.writelines("""    </body>\n""")

    htmlfile.writelines("""</html>\n""")

#    for line in lines:

#        htmlfile.writelines(line)
    htmlfile.close()
#    os.system("cat temphtml.txt > %s" % htmlfilename)
    return htmlfilename


def main():
    ##no need for fancy config file parsing, all files have same name in result directory
    outputdir=sys.argv[1]
    print 'Preparing html network/cytoscape output in directory '+outputdir
    html_prepare(outputdir)
    print 'Preparing result.html'
    result_html_prepare(outputdir)
    print 'Preparing summary html'
    summary_html(outputdir)
    print 'Preparing title html'
    title_html(outputdir)

if __name__=='__main__':
    main()
