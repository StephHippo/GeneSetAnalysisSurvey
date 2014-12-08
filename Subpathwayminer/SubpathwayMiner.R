# Parkinsons Data Set
parkinsonsgenes <- readLines(parkinsonsidpath)
# Annotate a set of genes to pathways
ann<-getAnn(parkinsonsgenes[1:10])
result <- printAnn(ann)
write.table(result,file="ParkinsonResultsPathways", col.names=NA,sep="\t")
# Annotate a set of genes to sub-pathways based on KEGG Orthology (KO) identifiers.
subGraphListKO<-getKcSubGraph(graphList=getDefaultKOUndirectedGraph())
annKO<-getKOAnn(parkinsonsgenes[1:10],graphList=subGraphListKO)
result <- printAnn(annKO)
write.table(result,file="ParkinsonResultsSubPathways", col.names=NA,sep="\t")

# Colorectal Cancer Data Set
colorectalgenes <- readLines(colorectalidpath)
ann<-getAnn(colorectalgenes[1:10])
result<-printAnn(ann)
write.table(result,file="ColorectalResultsPathways", col.names=NA,sep="\t")
subGraphListKO<-getKcSubGraph(graphList=getDefaultKOUndirectedGraph())
annKO<-getKOAnn(colorectalgenes[1:10],graphList=subGraphListKO)
result <- printAnn(annKO)
write.table(result,file="ColorectalResultsSubPathways", col.names=NA,sep="\t")

# Chromosome 13 Data Set
chr13genes=readLines(chromosome13path)
ann<-getAnn(chr13genes[1:4])
result <- printAnn(ann)
write.table(result,file="Chr13ResultsPathways", col.names=NA,sep="\t")
subGraphListKO<-getKcSubGraph(graphList=getDefaultKOUndirectedGraph())
annKO<-getKOAnn(chr13genes[1:4],graphList=subGraphListKO)
result <- printAnn(annKO)
write.table(result,file="Chr13ResultsSubPathways", col.names=NA,sep="\t")