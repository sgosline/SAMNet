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
    file.writelines("""<html><body><center><a href="http://fraenkel.mit.edu/samnetweb" target="_parent"><img src="../../../samnet/img/samnet_header.png"></a></center></body></html>""")

    file.close()

def summary_html(outdir,comcolors):
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
    <link rel="stylesheet" type="text/css" href="../../../main.css">

    </style>

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
<h1>   </h1>
<center><h1>Network Legend</h1><br>
   <table><tr><td colspan="2"><center><h3>Edge commodities</h3></center></td></tr>
        """%(tindex,nontindex,mrnatindex,mrnanontindex))
    ##add in colors here
    for comm in comcolors.keys():
        file.writelines('<tr><td><h4>'+comm+':</h4></td><td> <font color="'+comcolors[comm]+'">----------------</font></td></tr>\n')
    file.writelines("""
    <tr><td colspan="2"><center><h3>Node Types</h3></center></td></tr>
    <tr><td colspan="2">
        <img src="../../../samnet/img/legend.png" height="150" width="225"/>
</td></tr>
   </table></center>

<center>
<h1>  </h1>
<h1>Result Summary</h1>
<table>
<tr><td><h2>Input files</h2>
    <a href="../data/proteinWeights" >Protein weights inputs</a><br/>
    <a href="../data/exp" >mRNA expression inputs</a><br/>
    <a href="../data/tf2gene" >Transcription Factor to DNA interaction file</a>
    <br/><a href="../data/interactome" >Protein-Protein Interaction PKL file</a><br/>
</td></tr>
<tr><td><h2>Output Files</h2>
    <a href="samnetoutmultiComm_mcfs_symbol.sif">SAMNet network in SIF format (Cytoscape) </a><br/>
    <a href="samnetoutmultiComm_edge_commodity_symbol.eda">Edge Flow values </a><br/>
    <a href="samnetoutmultiComm_edge_type_symbol.eda">Edge Type values </a><br/>
    <a href="samnetoutmultiComm_node_flow_symbol.noa">Node Flow values </a><br/>
    <a href="samnetoutmultiComm_node_type_symbol.noa">Node Type values </a><br/>
    <a href="../samnet_CYTOEDA_Benjamini_1.0_sigDavidTerms.xls">DAVID terms for each commodity</a><br/></td></tr>
<tr><td> <h2>Input found in Solution</h2></td></tr>
<tr><td>
  <div id="visualization" style="width: 300px; height: 200px;"></div>
</td>
</tr>
<tr><td>
  <div id="visualizationdna" style="width: 300px; height: 200px;"></div>
</td></tr>
</table></center>
    </body></html>
""" )

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
#symbol='_symbol' ##set to '' if there are no symbols

#input_filename='samnetoutmultiComm'
def get_sif_filename(outdir,input_filename,symbol):
    return os.path.join(outdir,input_filename+'_mcfs'+symbol+'.sif')

def get_node_flow_noa_filename(outdir,input_filename,symbol):
    return os.path.join(outdir,input_filename+'_node_flow'+symbol+'.noa')

def get_node_type_noa_filename(outdir,input_filename,symbol):
    return os.path.join(outdir,input_filename+'_node_type'+symbol+'.noa')

def get_edge_commodity_eda_filename(outdir,input_filename,symbol):
    return os.path.join(outdir,input_filename+'_edge_commodity'+symbol+'.eda')

def get_edge_type_eda_filename(outdir,input_filename,symbol):
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
    dozencolors=["#003366","#339933","#999911","#FF9999","#9999FF","#66FFFF","#FF9933","#990000","#595959","#00FFFF","#CC9900","#996633"]
    commdict={}
    for ind,c in enumerate(commodities):
        commdict[c]=dozencolors[ind % len(dozencolors)]
    return commdict
    

def sifParser(outdir, input_filename,symbol):
    nodeset = set()
    edgedict=defaultdict(list)
    ##i dont think i need these three files
    pwfile = get_protweight_filename(outdir)
    exprfile = get_expr_filename(outdir)
    tf_file = get_tf_filename(outdir)
    
    ##this gets the symbol filename
    siffile = get_sif_filename(outdir,input_filename,symbol)
    
    ##but we really need the edges and nodes by commodity
    #edge_comm_filename=get_edge_commodity_eda_filename(outdir,input_filename)
    node_type_filename=get_node_type_noa_filename(outdir,input_filename,symbol)
    node_flow_filename=get_node_flow_noa_filename(outdir,input_filename,symbol)

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

'''
    <style>
    .title {
      width= 300;
      height=30;
    }
    .graph {
      width=250;
      height=300;
    }
    .legend {
      width=50;
      height=300;
    }
   </style>
'''
def iframe_result_html_prepare(outdir):
    '''
    Added new HTML prepare file that doesn't use frames since they have been depracated.  
    Fingers crossed!
    '''
    input_prefix='samnetoutmultiComm'
    input_filename=os.path.join(outdir,input_prefix)
    
    resultfilename = get_html_result(outdir, input_prefix)

    htmlfilename = get_html_filename(outdir,input_prefix)

    resultfile = open(resultfilename,'w')
    resultfile.writelines("""
    <html>
    <head>
    <title>SAMNetWeb - Results</title>
  </head>
  <body>
    <table>
    <tr><td colspan=2>
  <iframe class="title" src="title.html" height="350" width="100%" scrolling="no" seamless="yes"</iframe></td><tr>
  """)
    resultfile.writelines("""
    <tr><td><iframe class="graph" src="'+os.path.basename(htmlfilename)+'" scrolling="no" width="80%"></iframe></td>
    """)
    resultfile.writelines("""
    <td><iframe class="legend" src="CytoscapeWebControls.html" width="20%"></iframe></td></tr></table>
    </body>
    </html>
    """)
    resultfile.close()

def result_html_prepare(outdir):
    '''
    This is the main file, with the 3 inset frames
    '''
    input_prefix='samnetoutmultiComm'
    input_filename=os.path.join(outdir,input_prefix)
    
    resultfilename = get_html_result(outdir, input_prefix)

    htmlfilename = get_html_filename(outdir,input_prefix)

    resultfile = open(resultfilename,'w')
    resultfile.writelines("""
    <html>
    <title>SAMNetWeb - Results</title>
    <frameset rows="15%,85%">
    <frame src="title.html" name="north" scrolling="no" frameborder="0"/>
    <frameset cols="68%,32%">
    """)
    ##now add in html filename
    resultfile.write('<frame src="'+os.path.basename(htmlfilename)+'" frameborder="0" name="center"/>\n')
    ##now add in the rest
    resultfile.writelines("""
    <frame src="CytoscapeWebControls.html" frameborder="0" name="east"/>
    </frameset>
    
    </frameset>
    </html>

    """)
    resultfile.close()


def html_prepare(outdir):
    '''
    This is the cytoscape plugin page
    '''
    input_prefix='samnetoutmultiComm'
    input_filename=os.path.join(outdir,input_prefix)

    symbol='_symbol'
    if not os.path.exists(input_filename+'_mcfs'+symbol+'.sif'):
        symbol=''
    
    htmlfilename = get_html_filename(outdir,input_prefix)


    htmlfile = open(htmlfilename,'w')
    htmlfile.writelines("""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">\n""")
    htmlfile.writelines("""<html>\n""")

    htmlfile.writelines("""    <head>\n""")
    htmlfile.writelines("""        <title>SAMNet Web Output Page</title>\n""")
    htmlfile.writelines("""    <link rel="stylesheet" type="text/css" href="../../../main.css">\n""")   
    htmlfile.writelines("""    
      <style type="text/css">
      div.section {
        margin: 0 auto 20px auto;
        width: 60%;
        padding-left: 80px;
        border-top: 1px solid #cbcbcb;
        font-family: Verdana;
      }

      div.last {
        border-bottom: 1px solid #cbcbcb;
        padding-bottom: 20px;
      }

      textarea {
        vertical-align: middle;
      }

      div.options {
        margin-left: 40px;
      }
    </style>
    """)
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

    if not os.path.exists(input_filename+'_mcfs'+symbol+'.sif'):
        htmlfile.write("<h1>NO NETWORK CREATED, CHECK PARAMETERS AND TRY AGAIN</h1>")
        comcolors={}
    else:
        networkline,commodities = sifParser(outdir, input_prefix,symbol)
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
    return htmlfilename,comcolors


def main():
    ##no need for fancy config file parsing, all files have same name in result directory
    outputdir=sys.argv[1]
    print 'Preparing html network/cytoscape output in directory '+outputdir
    htmlfile,comcolors=html_prepare(outputdir)
    print 'Preparing result.html'
    result_html_prepare(outputdir)
    print 'Preparing summary html'
    summary_html(outputdir,comcolors)
    print 'Preparing title html'
    title_html(outputdir)

if __name__=='__main__':
    main()
