args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

VCFfile<-args[1]
outputFile<-args[2]
MyFunctions<-args[3]

source(MyFunctions)
cnVCF<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","ALL")

vcfHeader<-scan(VCFfile,sep = "\n",what = "character")
vcfHeader<-vcfHeader[grep("^##INFO",vcfHeader)]
cnInfo<-sapply(strsplit(vcfHeader,"=|,"),function(x) x[3])

vcf<-read.delim(VCFfile,header=FALSE,comment.char = "#",stringsAsFactors = FALSE)
colnames(vcf)<-cnVCF
resList<-sapply(strsplit(vcf$INFO,";"),strsplit,"=")
valInfo<-lapply(resList,function(x){sapply(x,function(x){ temp<-x[2]; names(temp)<-x[1]; temp})})
tabInfo<-data.frame(matrix(nrow = len(resList),ncol = len(cnInfo)))
colnames(tabInfo)<-cnInfo
vcf$INFO<-NULL
for(i in 1:length(valInfo)){tabInfo[i,names(valInfo[[i]])]<- valInfo[[i]]}
fastWrite(cbind(vcf,tabInfo),outputFile)
