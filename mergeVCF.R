args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

OUTPATH<-args[1]
MyFunctions<-args[2]

source(MyFunctions)

files<-list.files(paste0(OUTPATH,"/TSV"))
samplesName<- substr(files,1,nchar(files)-4)
tsvResList<-list()
for(i in 1:length(samplesName)){
	tsvRes<-fastRead(paste0(OUTPATH,"/TSV/",files[i]),row.names = NULL)
	tsvRes$sample<-samplesName[i]
	tsvResList[[samplesName[i]]]<-tsvRes
}
tsvRes<-do.call(rbind,tsvResList)

fastWrite(tsvRes,paste0(OUTPATH,"/results/merged.tsv"))

tsvRes<-tsvRes[tsvRes$DP>100,]
tsvRes$GENOTYPE<-substr(tsvRes$ALL,1,3)

require(ggplot2)
pdf(paste0(OUTPATH,"/results/barplotGenotype.pdf"),width = 16,height = 9)
ggplot(tsvRes,mapping = aes(x=sample,fill=GENOTYPE))+
	geom_bar()+theme_bw()+
	ylab("Number of variant with DP>100")+
	theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3))
dev.off()
	
