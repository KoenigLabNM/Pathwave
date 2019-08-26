#########################################################################################
#											#
#					PACKAGE USAGE:					#
#											#
#########################################################################################

#Mapping the expression data

#Read in gene-enzyme relations from KEGG's XML files in folder XML/ (this can take a while) 
keggXmlPath=paste(myPathWavePath,"XML/",sep="/")
keggXml=pwKEGGxml(keggXmlPath)


#Read the file agilentIDToEntrezGene.tab in folder Data/
agToEzFile=paste(myPathWavePath,"Data/agilentIDToEntrezGene.tab",sep="/")
agToEz=read.table(agToEzFile)
agToEz=as.matrix(agToEz)


#Read in the expression data, also in folder Data/
exprFile=paste(myPathWavePath,"Data/Neuroblastoma_vsn_expr_stage1_4amp.tab",sep="/")
x.org=read.table(exprFile,row.names=NULL)
x.org=as.matrix(x.org)


#Get Agilent-IDs:
agilent.probes=agToEz[,1]
entrez.IDs=as.numeric(agToEz[,2])

agilent.annotation=x.org[,1]
all.agilentIDs=unique(agilent.annotation)


#Read in the class information
clFile=paste(myPathWavePath,"Data/Neuroblastoma_vsn_expr_stage1_4amp.class",sep="/")
cl=read.table(clFile)


#Build class factor
y=NULL
names=colnames(x.org)[colnames(x.org)!="row.names"]
for(i in colnames(x.org)){
  y=c(y,as.vector(cl$V2)[i == cl$V1])
}
y=as.factor(y)


#Map gene expression values onto reactions - START
reac.entrez=NULL
for(maps in names(keggXml$genes)){
  toDel=NULL   

  #Check if more than one gene maps to a reaction
  duplicate.reactions=unique(names(keggXml$genes[[maps]])[duplicated(names(keggXml$genes[[maps]]))])
  for(i in duplicate.reactions){

    reac.entrez=rbind(reac.entrez,c(i,paste((keggXml$genes[[maps]])[grep(paste("^",i,"$",sep=""),names(keggXml$genes[[maps]]))],collapse="_")))    
    toDel=c(toDel,grep(i,names(keggXml$genes[[maps]])))
  } 
  if(is.null(toDel)){
    rest.reactions=names(keggXml$genes[[maps]])
  } else{
    rest.reactions=names(keggXml$genes[[maps]])[-toDel]
  }
  for(i in rest.reactions){
    reac.entrez=rbind(reac.entrez,c(i,(keggXml$genes[[maps]])[grep(paste("^",i,"$",sep=""),names(keggXml$genes[[maps]]))]))
  }
}
index=order(reac.entrez[,1])
reac.entrez=reac.entrez[index,]


#Map entrez IDs on reactions
x=NULL
x.names=NULL
for(i in 1:nrow(reac.entrez)){
  reac.entrezIDs=unlist(strsplit(reac.entrez[i,2],"_"))
  
  #Get agilent IDs
  reac.agilent=NULL
  for(j in reac.entrezIDs){
    help=agilent.probes[grep(paste("^",j,"$",sep=""),entrez.IDs,perl=TRUE)]
    
    if(length(help)>0){
      reac.agilent=c(reac.agilent,help)
    }
  }
  #Bind expression data together
  rows=agilent.annotation %in% reac.agilent
  
  if(length(rows[rows])>0){
    help=x.org[rows,2:ncol(x.org)]
    
    #Convert values into numbers (were read in as characters)
    #If more than one gene maps on a reaction: calculate mean value
    if(is.null(dim(help))){
      help=as.numeric(help)      
    } else{
      help=as.matrix(apply(help,2,as.numeric))
      help=apply(help,2,mean)
    }

    x=rbind(x,help)
    x.names=c(x.names,reac.entrez[i,1])
  }
}
rownames(x)=x.names
#Map gene expression values onto reactions - END


#Build adjacency matrices
adj=pwAdjMatrices(keggXml$reactions,keggXml$compounds)

#Get optimised grids
#Download the RData file optimalGridHSA.RData and load it
#load("optimalGridHSA.RData")

#or

#Get optimised grids from log files in folder OptLogFiles
lp.path=paste(myPathWavePath,"OptLogFiles/",sep="")
optimalGridHSA=pwOptGrids(lp.path,adj)



#Some reactions might be in x that are not in optimalGridHSA
all.reac=unlist(lapply(optimalGridHSA,function(entry){unique(as.vector(entry$M))[unique(as.vector(entry$M))!="0"]}))
all.reac=unique(all.reac)

toDel=which(is.na(match(rownames(x),all.reac)))
if(length(toDel)>0){
  x=x[-toDel,]
}


#Convert the expression values into Z-scores
x=t(apply(x,1,function(entry){(entry-mean(entry))/sd(entry)}))


print(paste("KEGG reactions:",length(all.reac),sep=" "))
print(paste("reactions with expression values:",nrow(x),sep=" "))


#Perform analysis
result=pathWave(x,y,optimalM=optimalGridHSA,pvalCutoff=0.01,genes=keggXml$genes)


#get table 1.
pval=NULL
score=NULL
all=NULL
down=NULL
up=NULL

for(i in names(result)){
  pval=c(pval,result[[i]]$p.value)
  score=c(score,result[[i]]$score)
  all=c(all,length(grep(i,rownames(x))))
  up=c(up,length(result[[i]]$reaction.regulation[result[[i]]$reaction.regulation<0]))
  down=c(down,length(result[[i]]$reaction.regulation[result[[i]]$reaction.regulation>0]))
}
res=cbind(up+down,all,up,down,pval,score)
rownames(res)=names(result)

print(res)