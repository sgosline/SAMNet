#!python
# by courtesy of HuangYi @ 20110424
##added by sgosline


def DAVIDenrich(listF, idType, bgF='', resF='', bgName = 'Background1',listName='List1', category = '', thd=0.1, ct=2):
    from suds.client import Client
    import os
    import logging
    logging.getLogger('suds.client').setLevel(logging.DEBUG)
   
    if len(listF) > 0 and os.path.exists(listF):
        inputListIds = ','.join(open(listF).read().split('\n'))
        print listName+' List loaded:'+inputListIds        
    else:
        print 'No list loaded.'
        raise

    flagBg = False
    if len(bgF) > 0 and os.path.exists(bgF):
        inputBgIds = ','.join(open(bgF).read().split('\n'))
        flagBg = True
        print 'Use file background.'
    else:
        print 'Use default background.'

    client = Client('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl')
    print 'User Authentication:',client.service.authenticate('sgosline@mit.edu')
    print idType
    listType = 0
    print 'Percentage mapped(list):', client.service.addList(inputListIds,idType,listName,listType)
    if flagBg:
        listType = 1
        print 'Percentage mapped(background):', client.service.addList(inputBgIds,idType,bgName,listType)

    print 'Use categories:', client.service.setCategories(category)
    chartReport = client.service.getChartReport(thd,ct)
    chartRow = len(chartReport)
    print 'Total chart records:',chartRow
    
    if len(resF) == 0 or not os.path.exists(resF):
        if flagBg:
            resF = listF + '.withBG.chartReport'
        else:
            resF = listF + '.chartReport'
    with open(resF, 'w') as fOut:
        fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
        for row in chartReport:
            rowDict = dict(row)
            categoryName = str(rowDict['categoryName'])
            termName='NA'
            termName = str(rowDict['termName'])
            listHits = str(rowDict['listHits'])
            percent = str(rowDict['percent'])
            ease = str(rowDict['ease'])
            Genes = str(rowDict['geneIds'])
            listTotals = str(rowDict['listTotals'])
            popHits = str(rowDict['popHits'])
            popTotals = str(rowDict['popTotals'])
            foldEnrichment = str(rowDict['foldEnrichment'])
            bonferroni = str(rowDict['bonferroni'])
            benjamini = str(rowDict['benjamini'])
            FDR = str(rowDict['afdr'])
            rowList = [categoryName,termName,listHits,percent,ease,Genes,listTotals,popHits,popTotals,foldEnrichment,bonferroni,benjamini,FDR]
            fOut.write('\t'.join(rowList)+'\n')
        print 'write file:', resF, 'finished!'

if __name__ == '__main__':
    from optparse import OptionParser
    parser=OptionParser()

    parser.add_option('--listName',default='name',dest='listName',help='Name of list')
    parser.add_option('--idType',default='OFFICIAL_GENE_SYMBOL',dest='idType',help='Registered ID type')
    opts,args=parser.parse_args()

    if len(args)==1:
        DAVIDenrich(listF = args[0], idType = opts.idType, listName = opts.listName, category = 'GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE')    

    elif len(args)==2:
        DAVIDenrich(listF = args[0], bgF=args[1],idType = opts.idType, listName = opts.listName, category = 'GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE')    
    elif len(args)==0:
        DAVIDenrich(listF = './list1.txt', idType = 'AFFYMETRIX_3PRIME_IVT_ID', listName = 'list1', category = 'GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,abcd,BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE')
 #   DAVIDenrich(listF = './list2.txt', idType = 'AFFYMETRIX_3PRIME_IVT_ID', listName = 'list2', category = 'abcd,BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE')
  #  DAVIDenrich(listF = './list2.txt', bgF = './list1.txt', idType = 'AFFYMETRIX_3PRIME_IVT_ID', bgName = 'Background1', listName = 'list2', category = 'abcd,BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE')

        
