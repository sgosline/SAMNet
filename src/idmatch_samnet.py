'''
idmatch_samnet.py
This file handles the file parsing to map various identifiers back and forth


Copyright (c) 2012 Sara JC Gosline
sgosline@mit.edu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

'''


import sys,os,csv,re,pickle
import types
from optparse import OptionParser

def main():
    '''
    handles various types of identifier matching
    '''
    parser=OptionParser()
    parser.add_option('--mappingFile',type='string',dest='map_file',default='',help='Mapping PKL file containing species-specific String identifier matching file. Defaults to human mapping.')
    parser.add_option('--aliasFile',type='string',dest='alias_file',default='',help='String protein alias file to be converted to pkl')
    parser.add_option("--outputdir",type='string',dest='output_dir',default='./',help='Directory to place pkl output file if aliasFile option is used')
    parser.add_option('--taxaId',type='string',dest='taxa',default='9606',help='Taxa id for species in question. DEFAULT is 9606 (human)')

    fpath=os.path.dirname(os.path.abspath( __file__ ))
    print fpath

    libpath=re.sub('src','lib',fpath)
    print libpath
    
    (opts,args)=parser.parse_args()
    allfiles=args

    if opts.map_file=='':
        if opts.taxa=='9606':
            geneDict=pickle.load(open(libpath+'/9606protein.aliases.v9.0_geneName.pkl','r'))
        elif opts.taxa=='10090':
            geneDict=pickle.load(open(libpath+'/10090protein.aliases.v9.0_geneName.pkl','r'))   
        else:
            geneDict=pickle.load(open(libpath+'/4932protein.aliases.v9.0_geneName.pkl','r'))   
    else:
        geneDict=pickle.load(open(opts.map_file,'r'))
    
    if len(allfiles)>0:
        for f in allfiles:
            parseSifFileFromStringToGeneName(open(f,'r'),os.path.dirname(f)+'/symbol_'+os.path.basename(f),geneDict)
    elif opts.alias_file!='':
        createStringToSymbolPkl(opts.alias_file,opts.taxa,opts.output_dir)

    

def parseSifFileFromStringToGeneName(sif_file,out_file,geneDict=''):
    '''
    reads an open sif file and replaces the STRING identifier with 
     a friendlier gene name
     '''
    lines=sif_file.readlines()
    new_file=csv.writer(open(out_file,'w'),delimiter='\t')
    for l in lines:
        arr=re.split('\t',l.strip())
        for a in range(len(arr)):

#            if '=' in arr[a]:##added to accomodate .noa and .eda files!!
#                pep=re.sub('10090.','',re.split('=',arr[a])[0])
#            else:
            pep=re.sub('4932.','',arr[a]) ##yeast names need to 
            pep=arr[a]
            arr[a]=pep.strip()
            if pep.strip() in geneDict.keys() and ('ENSP' in pep.strip() or 'ENSMUSP' in pep.strip() or len(re.sub('Y[A-Z][A-Z][0-9][0-9][0-9][W|C]','',pep.strip()))<1 or 'icrogid:' in pep.strip()):
                matches=geneDict[pep.strip()]
                if len(matches)>1 and len(re.sub('Y[A-Z][A-Z][0-9][0-9][0-9][W|C]','',pep.strip()))<1:#we're in yeast
                    for m in matches:
                        if re.sub('[A-Z][A-Z][A-Z][0-9]*','',m)=='':
                            break
                else:
                    m=matches[0]
                arr[a]=re.sub(pep.strip(),m,arr[a])
        new_file.writerow(arr)


def parseNodeAttrFileFromStringToGeneName(sif_file,out_file,geneDict='',isStringAtt=False):
    '''
    reads an open noa file and replaces with
     a friendlier gene name
     '''
    lines=sif_file.readlines()
    new_file=csv.writer(open(out_file,'w'),delimiter=' ')
    new_file.writerow(re.split(' ',lines[0].strip()))##stupid function needs an array
    attrdict={}

    for l in lines[1:]:
        arr=re.split(' = ',l.strip())
        pep=arr[0].strip()
        if pep in geneDict.keys() and ('ENSP' in pep or 'ENSMUSP'in pep or len(re.sub('Y[A-Z][A-Z][0-9][0-9][0-9][W|C]','',pep))<1 or 'icrogid:' in pep):  
            matches=geneDict[pep]
            if len(matches)>1 and len(re.sub('Y[A-Z][A-Z][0-9][0-9][0-9][W|C]','',pep))<1:#we're in yeast
                for m in matches:
                    if re.sub('[A-Z][A-Z][A-Z][0-9]*','',m)=='':
                        break
            else:
                m=matches[0]
            pep=m
        if pep not in attrdict.keys():
            attrdict[pep]=arr[1]
        elif isStringAtt:
            attrdict[pep]=attrdict[pep]+','+arr[1]
        else:
            print pep+' '+attrdict[pep]
            print pep+' '+arr[1]
            attrdict[pep]=str(float(attrdict[pep])+float(arr[1]))
            print attrdict[pep]
            
            
            
    for k in attrdict.keys():
        new_file.writerow([k,'=',attrdict[k]])


def parseAttrFileFromStringToGeneName(sif_file,out_file,geneDict=''):
    '''
    reads an open sif file and replaces the STRING identifier with 
     a friendlier gene name
     '''
    lines=sif_file.readlines()
    new_file=csv.writer(open(out_file,'w'),delimiter=' ')
    for l in lines:
        arr=re.split(' ',re.sub('\n','',l))
        for a in range(len(arr)):
            pep=arr[a]
            arr[a]=pep.strip()
            
            if pep.strip() in geneDict.keys() and ('ENSP' in pep.strip() or 'ENSMUSP'in pep.strip() or len(re.sub('Y[A-Z][A-Z][0-9][0-9][0-9][W|C]','',pep.strip()))<1 or 'icrogid:' in pep.strip()):  
                matches=geneDict[pep.strip()]
                if len(matches)>1 and len(re.sub('Y[A-Z][A-Z][0-9][0-9][0-9][W|C]','',pep.strip()))<1:#we're in yeast
                    for m in matches:
                        if re.sub('[A-Z][A-Z][A-Z][0-9]*','',m)=='':
                            break
                else:
                    m=matches[0]
                arr[a]=re.sub(pep.strip(),m,arr[a])

        new_file.writerow(arr)


def parseTabFileFromStringToGeneName(sif_file,out_file,geneDict=''):
    '''
    reads an open tab delimited file and replaces the STRING identifier with 
     a friendlier gene name
     '''
    lines=sif_file.readlines()
    if len(lines)==0:
        return(lines)
    new_file=csv.writer(open(out_file,'w'),delimiter='\t')
    new_file.writerow(re.split('\t',lines[0].strip()))
    new_data={}
    for l in lines[1:]:
        arr=re.split('\t',re.sub('\n','',l))
        pep=arr[0].strip()
        
        if pep in geneDict.keys() and ('ENSP' in pep.strip() or 'ENSMUSP'in pep.strip() or len(re.sub('Y[A-Z][A-Z][0-9][0-9][0-9][W|C]','',pep))<1 or 'icrogid:' in pep.strip()):  
            matches=geneDict[pep.strip()]
            if len(matches)>1 and len(re.sub('Y[A-Z][A-Z][0-9][0-9][0-9][W|C]','',pep.strip()))<1:#we're in yeast
                for m in matches:
                    if re.sub('[A-Z][A-Z][A-Z][0-9]*','',m)=='':
                        break
            else:
                m=matches[0]
            pep=m
        if pep not in new_data.keys():
            new_data[pep]=arr[1:]
        else: ##we have to resolve new identifiers
            new_arr=[]
            old_dat=new_data[pep]
            new_dat=arr[1:]
            #print old_dat
            #pKRTrint new_dat
            for i in range(len(new_dat)):
                new_arr.append(str(float(new_dat[i])+float(old_dat[i])))
            new_data[pep]=new_arr
    for k in new_data.keys():
        new_file.writerow([k]+new_data[k])

    


if __name__=='__main__':
    main()


