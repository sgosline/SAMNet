SAMNet

------------------------------------------------------------------------------------------------
Simultaneous Analysis of Multiple Networks (SAMNet) is a tool designed to identify key proteins undetected by high throughput biochemical experiments.  
------------------------------------------------------------------------------------------------
Contact: Sara JC Gosline sgosline@mit.edu
------------------------------------------------------------------------------------------------
Copyright (c) 2012-2016 Sara JC Gosline

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

------------------------------------------------------------------------------------------------
System Requirements:
------------------------------------------------------------------------------------------------
1- Python 2.5.2 or higher: http://python.org

2- Networkx 1.7 or higher: http://networkx.lanl.gov

3- NumPy: http://numpy.scipy.org

4- AMPL: http://www.ampl.com/DOWNLOADS/index.html

5- cplexamp or some other solver that will work with AMPL:
http://www.ampl.com/CPLEX/#academic

6-The Python SOAP library, suds: https://pypi.python.org/pypi/suds


------------------------------------------------------------------------------------------------
Getting started:
------------------------------------------------------------------------------------------------

1- Download code from source repository

2- Install ampl, cplexamp (or some other solver) and necessary python libraries

3- Investigate SAMNet run on the two datasets from Gosline et al. 2012.  The
yeast metal dataset from Jin et al. can be found in the ./data/yeast_metal
subdirectory.  The EMT dataset from Thomson et al. can be found in the
./data/human_emt data

4- Either run new results, or explore results included using the tools in the
src/cytoscape directory

5- Also investigate GO enrichment tools in the src/go_enrichment directory.

6- To run on your own data, you'll need to parse it into the appropriate file formats:
   -protein weight files, formated like the *.phen files in the sample data
   directories, 1 for each commodity
   -mRNA expressiond ata, formated like the *.txt files in the sample data
   directories
   -a protein-protein interaction network in the networkx format, or you can use
   one of the ones provided
   -a protein-DNA interaction network, in the format of the *.tfa files provided
   -for more details type python src/samnet.py --h



