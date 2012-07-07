###does some basic go enrichment tests
###This is an R script to handle go enrichment of various commodities compared to the whole network.
###Source this script to handle others
library(Category)
library(GOstats)
library(org.Sc.sgd.db)
library(KEGG.db)
library(org.Hs.eg.db)


##CHANGE THIS:
## Update list of commodities and species to reflect network
comnames=c('arsenic','cadmium','chromium','copper','mercury','silver','zinc','zeb1','fixed','tgfb','snail')
species.name='Yeast'
##
##

##general file management functions
removeSourceSinkNodes<-function(tab){
  remonodes=union(which(sapply(tab[,1],function(x) x%in%c('S','T',comnames,paste(comnames,'sink',sep='_')))),which(sapply(tab[,3],function(x) x%in%c('S','T',comnames,paste(comnames,'sink',sep='_')))))
  print(paste("removing ",length(remonodes),'edges'))
  return(tab[-remonodes,])
  
}

readFile<-function(filename){
  tab<-removeSourceSinkNodes(read.table(filename,sep='\t'))
  comm.prots=lapply(unique(tab[,2]),function(x) union(tab[which(tab[,2]==x),1],tab[which(tab[,2]==x),3]))
  names(comm.prots)<-unique(tab[,2])
  return(comm.prots)
}

##gets genes that are distinct to a single element of the list
getDistinctProts<-function(genelist){
  newlist<-sapply(names(genelist),function(g){
    ugenes=unique(genelist[[g]])
    for(a in setdiff(names(genelist),g)){
     # print(a)
      if(a%in%names(genelist))
        ugenes<-setdiff(ugenes,genelist[[a]])
#      print(ugenes)
    }
    print(paste('Got',length(ugenes),'unique genes for commodity',g,'out of',length(unique(genelist[[g]]))))
    return(ugenes)})

  newlist
}

getCommonNames<-function(gl,species=species.name){
  if(tolower(species)=='yeast'){
    map<-as.list(org.Sc.sgdGENENAME)
  }
  else{ ##assume its human
    map<-as.list(org.Hs.egSYMBOL2EG)###as.list(org.Hs.egSYMBOL)
#    print(head(map))
  }
    
  newgl <- sapply(gl,function(x,map){
    if(x%in%names(map)&&!is.na(map[[x]]))
      return(map[[x]])
    else
      return(x)},map)

 # print(head(gl))
 # print(head(newgl))
  return(newgl)

}

##computation of p-values here
computeHyperG<-function(foreground,allgenes=c(),ont='BP',pval=0.01,species=species.name){
  #first create a parameter object
  print(species)
  if(tolower(species)=='yeast'){
    if(length(allgenes)==0)
      allgenes=unique(unlist((as.list(org.Sc.sgdENSEMBL))))
    annote='org.Sc.sgd.db'
  }
  else{
    if(length(allgenes)==0)
      allgenes=unique(unlist((as.list(org.Hs.egGO2ALLEGS))))
    else
      allgenes<-getCommonNames(allgenes,species)
    foreground<-getCommonNames(foreground,species)
    annote='org.Hs.eg.db'
  }
  print(paste("All genes: ",paste(allgenes[1:10],collapse=',')))
  print(paste("Selected genes: ",paste(foreground[1:10],collapse=',')))
 # print(length(allgenes))
  pobj=new('GOHyperGParams',geneIds=foreground,universeGeneIds=allgenes,annotation=annote,ontology=ont,pvalueCutoff=pval,conditional=TRUE,testDirection='over')
  res<-hyperGTest(pobj)
  tab<-summary(res)##this contains all possible terms
 tab$fdr.pval=p.adjust(tab$Pvalue,"fdr")
  print(paste('Got',nrow(tab),'terms to correct'))

  restab<-tab[which(tab$Pval<pval),]
   print(head(restab))
  return(restab)
}

computeKeggHyperG<-function(foreground,allgenes=c(),pval=0.01,species=species.name){
  if(tolower(species)=='yeast'){
    annote='org.Sc.sgd.db'
    if(length(allgenes)==0)
      allgenes=unique(unlist((as.list(org.Sc.sgdGO2ALLORFS))))
   
  }
  else{
    annote='org.Hs.eg.db'
    if(length(allgenes)==0)
      allgenes=unique(unlist(as.list(org.Hs.egKEGG2ALLEGS)))
  }
  pobj=new('KEGGHyperGParams',geneIds=foreground,universeGeneIds=allgenes,annotation=annote,pvalueCutoff=1.0,testDirection='over')
  
  res<-hyperGTest(pobj)
  tab<-summary(res)##this contains all possible terms
  print(paste('Got',nrow(tab),'terms to correct'))
  tab$fdr.pval=p.adjust(tab$Pvalue,"fdr")
  restab<-tab[which(tab$Pval<pval),]
  return(restab)
###  return(summary(res))
}



## getUnCommonNames<-function(gl,species=species.name){
##   if(tolower(species)=='yeast')
##     map<-as.list(org.Sc.sgdGENENAME)
##   else ##assume its human
##     map<-as.list(org.Hs.egSYMBOL2EG)###as.list(org.Hs.egSYMBOL)
  
##   gl <- sapply(gl,function(x,map){
##     if(x%in%names(map)&&!is.na(map[[x]]))
##       return(map[[x]])
##     else
##       return(x)},map)
  
##   return(gl)

##}



###
#
# Default go table, computes GO terms for genes unique to a specific commodity
# using the full network as a background
#
##
buildGoTable<-function(siffile,filedescrip,go='BP',pval=0.01,type='GO',species=species.name){
  fullres<-NULL

  allpaths=readFile(siffile)
  protlist=getDistinctProts(allpaths)
  common=unique(unlist(allpaths))
  for(p in protlist)
    common=setdiff(common,p)
  print(paste("Got",length(common),'common genes'))
  allpaths$shared=common
  protlist$shared=common
    
  ##now write gene list to file
  genefile=c('Commodity','Unique Genes')
  commNames<-sapply(protlist,getCommonNames,species)
  names(commNames)<-names(protlist)
  for(a in names(protlist)){
#    pl<-getCommonNames(protlist[[a]])
    pl=commNames[[a]]
    genefile<-rbind(genefile,c(a,paste(pl,collapse=',')))
  }
#  genefile<-rbind(genefile,c("Shared",paste(getCommonNames(allpaths[["shared"]]),collapse=',')))
                                        #    fn=unlist(strsplit(g,split='/',fixed=T))
 # write.table(genefile,file=paste(filedescrip,'genes.txt',sep=''),row.names=FALSE,col.names=FALSE)
  res<-NULL
  if(tolower(species)=='human')
    protlist=commNames
  for(path in names(protlist)){
    if(length(protlist[[path]])==0)
      next
    if(type=='GO')
      tres<-computeHyperG(protlist[[path]],unique(unlist(protlist)),go,pval,species)
    else if(type=='KEGG')
      tres<-computeKeggHyperG(protlist[[path]],unique(unlist(protlist)),pval,species)
    print(paste('Got',nrow(tres),'terms for',path))
    if(nrow(tres)>0)
      res<-rbind(res,cbind(rep(path,nrow(tres)),tres))
  }
  if(!is.null(res)){
    print(head(res))
    fullres<-cbind(rep(filedescrip,nrow(res)),res)
  
    colnames(fullres)[1:2]<-c('siffile','commodity')
    filename=unlist(strsplit(siffile,split='/'))
    filename=gsub('.sif','',filename[length(filename)])
    for(cn in c('arsenic','cadmium','chromium','copper','mercury','silver','zinc')){
      filename=gsub(paste(cn,'_',sep=''),'',filename)
      filename=gsub(cn,'',filename)
    }
    fn=paste(type,'_TermsEnriched_distinctProts_netBG_At_p',pval,'_forNetwork_',filename,'.xls',sep='')
    write.table(fullres,fn,col.names=T,row.names=F,sep='\t')
#      ufullres=fullres[which(!duplicated(fullres$Term)),]
        dupterms=unique(fullres$Term[which(duplicated(fullres$Term))])
    uterms=setdiff(unique(fullres$Term),dupterms)
    ufullres=fullres[which(fullres$Term%in%uterms),]
    fn=paste(type,'_UNIQUE_TermsEnriched_distinctProts_netBG_At_p',pval,'_forNetwork_',filename,'.xls',sep='')
    if(nrow(ufullres)>0)
      write.table(ufullres,fn,col.names=T,row.names=F,sep='\t')
    return(fn)
  }
  else
    return('no terms')
}
##
#
# This function computes the enrichment for all terms in a specific commodity (including shared) 
# and computes their enrichment compared to the network BG
#
##
buildGoTableNonDistinct<-function(siffile,filedescrip,go='BP',pval=0.01,type='GO',species=species.name){
  fullres<-NULL
  
  allpaths=readFile(siffile)
  protlist=allpaths#getDistinctProts(allpaths)
  common=unique(unlist(allpaths))
 
  for(p in protlist)
    common=setdiff(common,p)
  print(paste("Got",length(common),'common genes'))
  allpaths$shared=common
  protlist$shared=common

  commNames<-sapply(protlist,getCommonNames,species)
  names(commNames)<-names(protlist)

  if(tolower(species)=='human')
    protlist=commNames
  
  res<-NULL
  for(path in names(protlist)){
    if(length(protlist[[path]])==0)
      next
    if(type=='GO')
      tres<-computeHyperG(protlist[[path]],unique(unlist(protlist)),go,pval,species)
    else if(type=='KEGG')
      tres<-computeKeggHyperG(protlist[[path]],unique(unlist(protlist)),pval,species)
    print(paste('Got',nrow(tres),'terms for',path))
    res<-rbind(res,cbind(rep(path,nrow(tres)),tres))
  }
  if(!is.null(res)){
    fullres<-rbind(fullres,cbind(rep(filedescrip,nrow(res)),res))
    
    colnames(fullres)[1:2]<-c('siffile','commodity')
    filename=unlist(strsplit(siffile,split='/'))
    filename=gsub('.sif','',filename[length(filename)])
    for(cn in c('arsenic','cadmium','chromium','copper','mercury','silver','zinc')){
      filename=gsub(paste(cn,'_',sep=''),'',filename)
      filename=gsub(cn,'',filename)
    }
 
    fn=paste(type,'_TermsEnrichedWithNetBG_for_nonDistinctGenes_At_p',pval,'_forNetwork_',filename,'.xls',sep='')
    write.table(fullres,fn,col.names=T,row.names=F,sep='\t')
#        ufullres=fullres[which(!duplicated(fullres$Term)),]
        dupterms=unique(fullres$Term[which(duplicated(fullres$Term))])
    uterms=setdiff(unique(fullres$Term),dupterms)
    ufullres=fullres[which(fullres$Term%in%uterms),]
    fn=paste(type,'_UNIQUE_TermsEnrichedWithNetBG_for_nonDistintGenes_At_p',pval,'_forNetwork_',filename,'.xls',sep='')
    if(nrow(ufullres)>0)
      write.table(ufullres,fn,col.names=T,row.names=F,sep='\t')
    return(fn)
  }
  else
    return('no terms')
}
##
#
# This function computes the enrichment for all terms in a specific commodity (including shared) 
# and computes their enrichment compared to the network BG
#
##
buildGoTableNonDistinctFullBG<-function(siffile,filedescrip,go='BP',pval=0.01,type='GO',species=species.name){
  fullres<-NULL
  
  allpaths=readFile(siffile)
  protlist=allpaths#getDistinctProts(allpaths)
  common=unique(unlist(allpaths))
 
  for(p in protlist)
    common=setdiff(common,p)
  print(paste("Got",length(common),'common genes'))
  allpaths$shared=common
  protlist$shared=common

  commNames<-sapply(protlist,getCommonNames,species)
  names(commNames)<-names(protlist)

  if(tolower(species)=='human')
    protlist=commNames
  
  if(tolower(species)=='yeast'){
    univ=unique(unlist((as.list(org.Sc.sgdGO2ALLORFS))))
  }
  else{
    if(type=='GO'){
      univ=unique(unlist((as.list(org.Hs.egGO2ALLEGS))))
    }
    else if(type=='KEGG'){
      univ=unique(unlist(as.list(org.Hs.egKEGG2ALLEGS)))##?? check this      
     # tres<-computeHyperG(entrezlist[[path]],univ,go,pval)
    }

    #tres<-computeKeggHyperG(entrezlist[[path]],univ,pval)
  }
  res<-NULL
  for(path in names(protlist)){
    if(length(protlist[[path]])==0)
      next
    if(type=='GO')
      tres<-computeHyperG(protlist[[path]],univ,go,pval,species)
    else if(type=='KEGG')
      tres<-computeKeggHyperG(protlist[[path]],univ,pval,species)
    print(paste('Got',nrow(tres),'terms for',path))
    res<-rbind(res,cbind(rep(path,nrow(tres)),tres))
  }
  if(!is.null(res)){
    fullres<-rbind(fullres,cbind(rep(filedescrip,nrow(res)),res))
    
    colnames(fullres)[1:2]<-c('siffile','commodity')
    filename=unlist(strsplit(siffile,split='/'))
    filename=gsub('.sif','',filename[length(filename)])
        for(cn in c('arsenic','cadmium','chromium','copper','mercury','silver','zinc')){
      filename=gsub(paste(cn,'_',sep=''),'',filename)
      filename=gsub(cn,'',filename)
    }
 
    fn=paste(type,'_TermsEnrichedWithFullBG_for_nonDistinctGenes_At_p',pval,'_forNetwork_',filename,'.xls',sep='')
    write.table(fullres,fn,col.names=T,row.names=F,sep='\t')
    dupterms=unique(fullres$Term[which(duplicated(fullres$Term))])
    uterms=setdiff(unique(fullres$Term),dupterms)
    ufullres=fullres[which(fullres$Term%in%uterms),]
      
    fn=paste(type,'_UNIQUE_TermsEnrichedWithFullBGfor_nonDistinctGenes_at_p',pval,'_forNetwork_',filename,'.xls',sep='')
    if(nrow(ufullres)>0)
      write.table(ufullres,fn,col.names=T,row.names=F,sep='\t')
    return(fn)
  }
  else
    return('no terms')
}
#
#
#This function collects proteins unique to each commodity and computes go/kegg enrichment for each
#compared to the full BG for that species
#

buildGoTableWithFullBG<-function(siffile,filedescrip,go='BP',type='GO',pval=0.01,species=species.name){
  fullres<-NULL
  
  allpaths=readFile(siffile)
  protlist=getDistinctProts(allpaths)
  
  common=unique(unlist(allpaths))
  for(p in protlist)
    common=setdiff(common,p)
  allpaths$shared=common
  protlist$shared=common

  commNames<-sapply(protlist,getCommonNames)
  names(commNames)<-names(protlist)

  if(tolower(species)=='human')
    protlist=commNames

  
  res<-NULL
  if(tolower(species)=='yeast'){
      univ=unique(unlist((as.list(org.Sc.sgdGO2ALLORFS))))
  }
  else{
    if(type=='GO'){
      univ=unique(unlist((as.list(org.Hs.egGO2ALLEGS))))
    }
    else if(type=='KEGG'){
      univ=unique(unlist(as.list(org.Hs.egKEGG2ALLEGS)))##?? check this      
     # tres<-computeHyperG(entrezlist[[path]],univ,go,pval)
    }

    #tres<-computeKeggHyperG(entrezlist[[path]],univ,pval)
  }
  for(path in names(protlist)){
    if(length(protlist[[path]])==0)
      next
    if(type=='GO')
      tres<-computeHyperG(protlist[[path]],univ,go,pval,species)
    else if(type=='KEGG')
      tres<-computeKeggHyperG(protlist[[path]],univ,pval,species)
    
    print(paste('Got',nrow(tres),'terms for',path))
    res<-rbind(res,cbind(rep(path,nrow(tres)),tres))
  }
  if(!is.null(res)){
    fullres<-rbind(fullres,cbind(rep(filedescrip,nrow(res)),res))
    
    colnames(fullres)[1:2]<-c('siffile','commodity')
    filename=unlist(strsplit(siffile,split='/'))
    filename=gsub('.sif','',filename[length(filename)])
        for(cn in c('arsenic','cadmium','chromium','copper','mercury','silver','zinc')){
      filename=gsub(paste(cn,'_',sep=''),'',filename)
      filename=gsub(cn,'',filename)
    }

    
    fn=paste(type,'_TermsEnrichedWithFullBGAt_p',pval,'_forNetwork_',filename,'.xls',sep='')
    write.table(fullres,fn,col.names=T,row.names=F,sep='\t')
    
    
    dupterms=unique(fullres$Term[which(duplicated(fullres$Term))])
    uterms=setdiff(unique(fullres$Term),dupterms)
    ufullres=fullres[which(fullres$Term%in%uterms),]

    fn=paste(type,'_UNIQUE_TermsEnrichedWithFullBGAt_p',pval,'_forNetwork_',filename,'.xls',sep='')
    if(nrow(ufullres)>0)
      write.table(ufullres,fn,col.names=T,row.names=F,sep='\t')
    
    return(fn)
  }
  else
    return("no terms")
  
}


