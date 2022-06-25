#!/usr/bin/R


lire<-function(x, sep = "\t",header=TRUE,row.names = 1,quote="",stringsAsFactors=TRUE,as.matrix=FALSE,...){
	res<-read.table(file=x,sep=sep,header=header,row.names = row.names,quote=quote,stringsAsFactors=stringsAsFactors,...)
  if(as.matrix){
    as.matrix(res)
  }else{
    res
  }
}


ecrire<-function(x,file="default.tsv",headRow="Name",row.names=TRUE,col.names=TRUE, dec=".",...){
	options(warn=-1) #Supress unecessary warning about append 
	if(row.names && col.names){
		write.table(x = paste0(headRow,"\t"),file = file,sep = "\t",eol="",quote=F,row.names=F,col.names=F)
		write.table(x=x,file=file,sep="\t", row.names = T, col.names = T, quote = FALSE,append=T,dec=dec,...)
	}else{
		write.table(x=x,file=file,sep="\t", row.names = row.names, col.names = col.names, quote = FALSE, dec=dec,...)
	}
	options(warn=0)
}


fastRead <- function(fileName, sep = '\t',row.names = 1,as.matrix=FALSE,stringsAsFactors=FALSE,...){
	require(data.table)
	dat <- as.data.frame(data.table::fread(fileName,stringsAsFactors=stringsAsFactors, sep = sep,...))
	if(!is.null(row.names)){
	  rownames(dat) <- dat[,row.names]
	  dat <- dat[,-row.names,drop=FALSE]
	}
	if(as.matrix) dat<-as.matrix(dat)
	return(dat)
}

fastWrite <- function(x, fileName = "default.tsv", headRow="Name",row.names=TRUE,col.names=TRUE, dec=".",sep="\t",...) {
	require(data.table)
	if(is.null(rownames(x))) row.names<-FALSE
	if(is.null(colnames(x))) col.names<-FALSE
	
	if(row.names){
		x=cbind(rownames(x),x)
		colnames(x)[1]<-headRow
	}
	fwrite(x=data.frame(x),file=fileName,sep=sep, row.names = FALSE, col.names = col.names, quote = FALSE, dec=dec,...)
}


write.vectorList <- function(list, filename,sep="\t",list.names=TRUE,vector.names=TRUE) {
	if( (!is.list(list))) stop("list must be a list")
	sink(filename)
	for(i in seq_along(list)){
		if(! (is.null(names(list)) | (!list.names))) cat(names(list[i]),"\n",sep = "")
		element<-list[[i]]
		if(! (is.vector(element) | is.factor(element))) stop("each element of the list should be a vector")
		if(! (is.null(names(element)) | (!vector.names))) cat(paste0(names(element),collapse = sep),"\n",sep = "")
		cat(paste0(as.character(element),collapse = sep),"\n",sep = "")
	}
	sink()
}

read.vectorList<-function(file,sep="\t"){
	require(stringr)
	con<-file(file)
	txt<-readLines(con)
	elNames<-str_remove_all(txt[seq(1,length(txt),2)],"\t")
	res<-strsplit(txt[seq(2,length(txt),2)],"\t")
	names(res)<-elNames
	close(con)
	res
}

get_rand_state <- function() {
  # Using `get0()` here to have `NULL` output in case object doesn't exist.
  # Also using `inherits = FALSE` to get value exactly from global environment
  # and not from one of its parent.
  get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
}

set_rand_state <- function(state) {
  # Assigning `NULL` state might lead to unwanted consequences
  if (!is.null(state)) {
    assign(".Random.seed", state, envir = .GlobalEnv, inherits = FALSE)
  }
}

gmean<-function(x, keepZero=FALSE){ #geometrical mean
	if(sum(x)==0) return(0)
	if(!keepZero){
		x<-x[x!=0]
	}else{
		if(length(which(x==0))>0) return(0)
	}
	return( exp( sum(log(x))/length(x) ) )
}


plotText<-function(cols,label=colnames(cols),cex=.4,xlab="",ylab="",main=""){
	plot(cols,type="n",xlab=xlab,ylab=ylab,main=main)
	text(cols,labels=label,cex=cex)
}

plotMatrixText<-function(mat,main="",cex=1,rowtitle="",rownames=TRUE,colnames=TRUE,insertBlank=c()){
	#Converts a character matrix into a plot on the graphic output
	matTxt<-matrix("",nrow(mat),ncol(mat))
	for(j in 1:ncol(mat)) matTxt[,j]<-as.character(mat[,j])
	for(j in insertBlank) matTxt<-cbind(matTxt[,1:j],rep("",nrow(matTxt)),matTxt[,(j+1):ncol(matTxt)])
	if(colnames){
		cnMat<-colnames(mat)
		for(j in insertBlank) cnMat<-c(cnMat[1:j],"",cnMat[(j+1):length(cnMat)])
		matTxt<-rbind(cnMat,matTxt)
	}
	if(rownames) matTxt<-cbind(c(rowtitle,rownames(mat)),matTxt)
	#coord calculation
	coordMatrixX<- matrix(0,nrow(matTxt),ncol(matTxt))
	coordMatrixY<-coordMatrixX
	for(i in nrow(coordMatrixX):1){
	  for(j in ncol(coordMatrixY):1){
		coordMatrixX[i,j]<-nrow(coordMatrixX)-i
		coordMatrixY[i,j]<-j
	  }
	}
	plot(1,xlim=c(1,ncol(matTxt)+2),ylim=c(0,nrow(matTxt)),type="n", axes = F, main=main,ylab="",xlab="")
	text(coordMatrixY,coordMatrixX,labels = matTxt,cex=cex,adj = c(0,0))
}

genColWithGrep<-function(nameVector,patternVector){
	#give a color vector from a character vector with a pattern vector
	cols=rainbow(length(patternVector))
	colVector<-rep("#FFFFFFFF",length(nameVector))
	i<-1
	for(pattern in patternVector){
		colVector[grep(pattern,nameVector)]<-cols[i]
		i<-i+1
	}
	return(colVector)
}

genColWithFactors<-function(factorVector){
	levels<-levels(as.factor(factorVector))
	cols=rainbow(length(levels))
	names(cols)<-levels
	return(cols)
	
}

genCOlWithTab<-function(x){
  col<-list()
  for(l in names(x)) col[[l]]<-rgb(x[[l]]$r/255,x[[l]]$g/255,x[[l]]$b/255,names = rownames(x[[l]]))
  return(col)
}

linearColorScale = function(x,colMin="white",colMax="red"){
	require(circlize)
	colorRamp2(c(min(x),max(x)),colors = c(colMin,colMax))(x)
}


cv<-function(x){ #Coefficient of variation
	return(sd(x)/mean(x));
}

cv2<-function(x){ #Coefficient of variation of Lanner
	return(sd(x)^2/mean(x)^2);
}

se<-function(x){ #Standard mean error
	return(sd(x)/sqrt(length(x)));
}

uncenter<-function(x){
	#transform vector to have no negative value
	return(x+abs(min(x)));
}

takefirst<-function(x,returnIndex=FALSE){
	uniqDat<-unique(x)
	caseUniq<-c()
	for(i in uniqDat) caseUniq<-c(caseUniq,which(i==x)[1])
	if(returnIndex){
		return(caseUniq)
	}else{
		return(x[caseUniq])
	}
}


LevenshteinDist<-function(vectA,vectB){
	lenA<-length(vectA);
	ncols<-lenA+1;
	lenB<-length(vectB);
	nrows<-lenB+1;
	CompMat<-data.frame(matrix(0,nrows,ncols));
	print(dim(CompMat))
	CompMat[,1]<-0:lenB;
	CompMat[1,]<-0:lenA;
	print("mat initialized")
	for(i in 2:nrows){
		for(j in 2:ncols){
			cost<-if(vectA[j-1]==vectB[i-1]) 0 else 1;
			CompMat[i,j]<-min(CompMat[i-1,j]+1, CompMat[i,j-1]+1, CompMat[i-1,j-1]+cost)
		}
		print("i=")
		print(i)
	}
	return(CompMat[nrows,ncols])
}

#Aggregation de ligne/colonne selon une variable qualitative
#obsolete, see aggregate
aggregCols<-function(dataframe,vector,fun=mean){ 
  vector<-as.factor(as.character(vector))
  #pour forcer le mise à jour de l'attribut levels
  lvl<-levels(vector)
  res<-matrix(0,nrow(dataframe),length(lvl))
  rownames(res)<-rownames(dataframe)
  colnames(res)<-lvl
  for(i in lvl){
	if(length(which(vector==i))>1){
		res[,i]<-apply(dataframe[,which(vector==i)],1,fun)
	}
	else{
		res[,i]<-dataframe[,which(vector==i)]
	}
  }
  return(res)
}

aggregRows<-function(dataframe,vector,fun=mean){
  vector<-as.factor(as.character(vector))
  lvl<-levels(vector)
  res<-matrix(0,length(lvl),ncol(dataframe))
  colnames(res)<-colnames(dataframe)
  rownames(res)<-lvl
  for(i in lvl){
	if(length(which(vector==i))>1){
		res[i,]<-apply(dataframe[which(vector==i),,drop=FALSE],2,fun)
	}else if(length(which(vector==i))==1){
		res[i,]<-unlist(dataframe[which(vector==i),,drop=FALSE])
	}
  }
  return(res)
}

autoGparFontSizeMatrix<-function(n,...){ #Calcule automatiquement la taille de police selon le nombre de colonnes ou lignes (empirique)
	n=max(n,50)
	n=min(n,1000)
	return(gpar(fontsize=1/n*600,...))
}

ConvertKey<-function(keyList, tabKey,colOldKey=1,colNewKey=2){
	#@param keyList: vector of id to be converted
	#@param tabKey : dataframe of x columns with at least two contain old and new keys correspondance
	#@param colOldKey : index of old keys column in tabKey 
	#@param colNewKey : index of new keys column in tabKey 
	hashCorr<-tabKey[,colNewKey]
	names(hashCorr)<-tabKey[,colOldKey]
	returned<-hashCorr[keyList]
	names(returned)<-NULL
	return(as.character(returned))
}

ConvertKeyMatrix<-function(tab,tabKey,colOldKey=1,colNewKey=2,first=TRUE,dim=1,fun){
	if(dim==2) tab<-t(tab)
	keyList<-rownames(tab)
	newkey<-ConvertKey(keyList, tabKey,colOldKey,colNewKey)
	names(newkey)<-as.character(1:length(newkey))
	if(first){
	  newkey<-newkey[which(!is.na(newkey))]
	  newkey<-takefirst(newkey)
	  tab<-tab[as.numeric(names(newkey)),]
	  rownames(tab)<-newkey
	}else{
	  tab<-tab[which(!is.na(newkey)),]
	  newkey<-newkey[which(!is.na(newkey))]
	  tab<-aggregRows(tab,newkey,fun=fun)
	}
	if(dim==2) tab<-t(tab)
	return(tab)
}

supprNAnames<-function(x ,side=1){
	#Vector case
	if(is.null(dim(x))){
		return(x[which(!is.na(names(x)))]);
	}
	#Col case
	if(side==2){
		return(x[,which(!is.na(colnames(x)))]);
	}
	#Row case
	else{
		return(x[which(!is.na(rownames(x))),]);
	}
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

mobileMean<-function(datVector, step=2){
#Lissage de signal à partir de moyenne mobile
  pond<-datVector
  for(i in 1:length(datVector)){
    pond[i]=mean(datVector[max(1,(i-step)):min((i+step),length(datVector))])
  }
  return(pond)
}

OrderBestUpDown<-function(x){ #-3,4,-2,2,6 --> 6,4,2,-3,-2
#x : named vector with positive and negative values  
	up<-sort(x[which(x>0)],decreasing =TRUE)
	down<-sort(x[which(x<=0)])
	return(c(up,down))
}


rowScale<-function(data,center=TRUE,scaled=FALSE){
	data<-t(data)
	data<-t(scale(data,center=center,scale=scaled))
	return(data)
}

corrDist<-function(x,method="pearson") return(as.dist((1 - cor(Matrix::t(x),method=method))/2))
corrDistBicor<-function(x){
  require("WGCNA")
  resDist<-as.dist((1 - suppressWarnings(bicor(Matrix::t(x))))/2)
}

chead<-function(x, n=5){
	print(x[1:n,1:n])
}

inflaslope<-function(x, ordering=TRUE, decreasing = FALSE){
#Calculer point d'ébouli des valeur, marche pas
	if(ordering) x<-sort(x, decreasing=decreasing)
	l1<-length(x)
	t1<-cumsum(x)
	t2<-abs(t1[1:(l1-1)]-t1[2:l1])[1:(l1-1)]
	t3<-which(t2>sd(x))[1]
	return(t3)
}

majR<-function(){
	# installing/loading the package:
	if(!require(installr)) {
	install.packages("installr"); require(installr)} #load / install+load installr
	# using the package:
	updateR()
}

# @definition Logical combination af a vector with custom operator #Maj : use sum !!!
combineLogical<-function(x,operator="&"){
	if(!is.vector(x))  stop("error, x must be a vector")
	if(!is.logical(x)) stop("error, x must be logical")
	if(!operator%in%c("&","|")) stop("error, operator must be equal to '&' or '|'")
	res<-sum(x)
	lg<-length(x)
	if(res==0) return(FALSE)
	if(operator=="&" & res<lg) return(FALSE)
	return(TRUE)
}

Mode <- function(x) {
# code from https://github.com/benmarwick/LaplacesDemon/blob/master/R/Mode.R
### Initial Checks
	if(missing(x)) stop("The x argument is required.")
	if(!is.vector(x)) x <- as.vector(x)
	x <- x[is.finite(x)]
	### Discrete
	if(all(x == round(x))) {
		Mode <- as.numeric(names(which.max(table(x))))
	} else {
		### Continuous (using kernel density)
		x <- as.vector(as.numeric(as.character(x)))
		kde <- density(x)
		Mode <- kde$x[kde$y == max(kde$y)]}
	return(Mode)
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

rmNA<- function(x) x<-x[!is.na(x)]

convertColorAdd2Sub<-function(red,green,blue){
	newRed=rowMeans(cbind(1-green,1-blue))
	newGreen=rowMeans(cbind(1-red,1-blue))
	newBlue=rowMeans(cbind(1-red,1-green))
	return(data.frame("red"=newRed,"green"=newGreen,"blue"=newBlue))
}

#create empty dataframe
emptyDF<-function(row.names,col.names,defaultValue=NA){
	tempMat<-matrix(data = defaultValue,nrow=length(row.names),ncol = length(col.names),dimnames = list(row.names,col.names))
	data.frame(tempMat)
}

#copy paste ready vector
copyReadyVector<-function(x){
	paste0("c('",paste0(x,collapse = "','"),"')")
}

#Better make.unique
make.unique2<-function(sample.name,sep="."){
	counts=table(as.factor(sample.name))
	nmRep<-sapply(as.list(counts),function(x) 1:x)
	paste0(rep(names(nmRep),counts),sep,as.character(unlist(nmRep,use.names = F)))[order(sample.name)[order(sample.name)]]
}

#str split with choosen element
strsplitNth<-function(x, split, n=1, fixed=FALSE, perl=FALSE, useBytes=FALSE){
	res<-strsplit(x, split, fixed, perl, useBytes)
	sapply(res,function(el){ el[n] })
}

matrixCoord1D_2D<-function(x,Mat){
	matDim<-dim(Mat)
	c((x-1) %% matDim[1] +1 , (x-1) %/% matDim[1] + 1)
}

grey2red = function(n, base, gamma)
{
  red = seq(from=base^gamma, to=1, length.out = n)^(1/gamma)
  green = blue = seq(from = base^gamma, to=0, length.out = n)^(1/gamma);
  col = rgb(red, green, blue, maxColorValue = 1); 
}

formatNumber2Character<-function(x,digit=max(nchar(as.character(x)))){
	x<-as.character(x)
	sapply(as.list(x),function(el){ paste0(paste0(rep("0",2-nchar(el)),collapse = ""),el) })
}

matrixFromDimnames<-function(row,col,value=0){
	matrix(value,ncol=length(col),nrow=length(row),dimnames=list(row,col))
}

pickChangeInVectorOrder<-function(x){
	c(x[1],x[2:length(x)] [x[1:(length(x)-1)]!=x[2:length(x)]])
	#suppressWarnings(x[x!=x[2:length(x)]]) method 2
}

mostDistantColor<-function(n,rlimits=c(0,1),glimits=c(0,1),blimits=c(0,1),resolution=50,seed=666){
	old_state <- get_rand_state()
  	on.exit(set_rand_state(old_state))
	
	set.seed(seed)
	colorSpace<-expand.grid(
		seq(rlimits[1],rlimits[2],length.out = resolution),
		seq(glimits[1],glimits[2],length.out = resolution),
		seq(blimits[1],blimits[2],length.out = resolution)
	)
	suppressWarnings(kcenters<-kmeans(x = colorSpace,centers = n)$centers)
	as.vector(apply(kcenters,1,function(x){
		rgb(red = x[1],green = x[2],blue  =x[3])
	}))
}

mostDistantColor2<-function(n, colorspace = "rainbow", cvd = c("protan", "deutan","tritan"), cvd_severity = 0, n_threads = NULL){
	require(qualpalr)
	if(n==1) return("#000000")
	qualpal(n=n, colorspace = colorspace, cvd = cvd, cvd_severity = cvd_severity, n_threads = n_threads)$hex
}

alphabetOrderByAdding0<-function(x){
	maxChar<-max(nchar(x))
	zeros<-sapply(x,function(el){
		paste0(el<-rep("0",maxChar-nchar(el)),collapse="")
	})
	paste0(zeros,x)
}

factorAsStrings<-function(dt){
	for(i in 1:ncol(dt)){
		if(is.factor(dt[,i])) dt[,i]<-as.character(dt[,i])
	}
	dt
}

factorToVectorList<-function(factorValues,factorNames=NULL){
	if(is.null(factorNames)) factorNames<-names(factorValues)
	res<-lapply(levels(factorValues),function(x) factorNames[factorValues==x])
	names(res)<-levels(factorValues)
	res
}

VectorListToFactor<-function(listOfVector){
	res<-factor(unlist(lapply(seq_along(listOfVector),function(i) rep(names(listOfVector)[i],length(listOfVector[[i]])))),
								 levels=names(listOfVector))
	names(res)<-unlist(listOfVector)
	res
}

scaleRangePerRow<-function(correctedExpr,oldExpr){ # minimum of gene expression is same before batch effect correction
	correctedExpr<-as.matrix(correctedExpr)
	oldExpr<-as.matrix(oldExpr)
	correctedExpr-apply(correctedExpr,1,min)+apply(oldExpr,1,min)
}

evaluateThresholdMetric<-function(metric,nmadArgumentName,nmads=NULL,nmads.upper=NULL,
																	nmads.lower=NULL,lower.limit=NULL,upper.limit=NULL){
	if(is.null(lower.limit) | is.null(upper.limit)){
		med <-median(metric, na.rm = TRUE)
		mad <-mad(metric, center = med, na.rm = TRUE)
		if(is.null(lower.limit)){
			if(is.null(nmads.lower)){
				if(is.null(nmads)) stop(paste0("You must give at least the 'nmads' argument."))
				nmads.lower=nmads
			}
			lower.limit=med-(nmads.lower * mad)
		}
		if(is.null(upper.limit)){
			if(is.null(nmads.upper)){
				if(is.null(nmads)) stop(paste0("You must give at least the 'nmads' argument."))
				nmads.upper=nmads
			}
			upper.limit=med+(nmads.upper * mad)
		}
	}
	c(lower.limit,upper.limit)
}

linearScale <- function(vals,newRange,returnFunction = TRUE) {
	if(!is.numeric(vals)) stop("x should be a vector of numerics")
	if(length(newRange)!=2 | !is.numeric(newRange)) stop("newRange should be a vector of 2 numerics")
	
	oldMin<-min(vals)
	oldMax<-max(vals)
	newMin<-newRange[1]
	newMax<-newRange[2]
	
	mfac<-(newMax-newMin)/(oldMax-oldMin)
	scaleFun<-function(x) newMin+(x-oldMin)*mfac

	if(returnFunction){
		scaleFun
	}else{
		scaleFun(vals)
	}
}

matDist<-function(M1,M2,method="euclidean"){
	if(!is.matrix(M1)) M1<-as.matrix(M1)
	if(!is.matrix(M2)) M2<-as.matrix(M2)
	
	res<-matrix(apply(expand.grid(1:nrow(M1),1:nrow(M2)),1,function(row){
		as.vector(dist(rbind(M1[row[1],],M2[row[2],]),method = method))
	}),nrow =nrow(M1))
	
	rownames(res)<-rownames(M1)
	colnames(res)<-rownames(M2)
	
	res
}

matDistEuclidean<-function(M1,M2){
	if(!is.matrix(M1)) M1<-as.matrix(M1)
	if(!is.matrix(M2)) M2<-as.matrix(M2)
	
	res<-sqrt(do.call("+",lapply(1:ncol(M1),function(i){
		outer(M1[,i],M2[,i],"-")^2
	})))
	
	rownames(res)<-rownames(M1)
	colnames(res)<-rownames(M2)
	
	res
}

pointdensity.nrd<-function(mat,eps=1){
	if(!is.matrix(mat)) mat<-as.matrix(mat)
	require(dbscan)
	require(MASS)
	pointdensity(apply(mat,2,function(d) d/bandwidth.nrd(d)),eps = eps,type = "density")
}


add_modify_aes <- function(mapping, ...) {
	ggplot2:::rename_aes(modifyList(mapping, ...))  
}

log10plus1<-function(){
	require(scales)
	scales::trans_new(name = "log10plus1",transform = function(x) log10(x+1),inverse = function(x) 10^x-1 ,domain=c(0,Inf))
}

intersectionEnrichment<-function(isInGroupMat){
	universe<-nrow(isInGroupMat)
	expected<-prod(apply(isInGroupMat,2,sum)/universe)*universe
	real<-sum(apply(isInGroupMat,1,function(x) sum(x))==ncol(isInGroupMat))
	real/expected
}

formatAnnotFromMeta<-function(annotDataFrame, metaAnnot){
	colorScales=list()
	for(feature in rn(metaAnnot)){
		if(! feature %in% colnames(annotDataFrame))
		annotDataFrame[,feature]<- do.call(paste0("as.",metaAnnot[feature,"Type"]),list(x=annotDataFrame[,feature]))
		if(metaAnnot[feature,"colorScale"] != ""){
			colorScale=strsplit(str_remove_all(metaAnnot[feature,"colorScale"]," "),split = ",")[[1]]
			splitted=sapply(colorScale,strsplit,"=")
			colorScales[[feature]]<-sapply(splitted,function(x) x[2])
			names(colorScales[[feature]])<-sapply(splitted,function(x) x[1])
			if(metaAnnot[feature,"Type"]=="factor"){
				annotDataFrame[,feature]<-factor(annotDataFrame[,feature],levels = names(colorScales[[feature]]))
			}
		}
	}
	attr(annotDataFrame,"colorScales")<-colorScales
	annotDataFrame
}

greps<-function(pattern, x, ignore.case = FALSE, perl = FALSE, value = FALSE,
								fixed = FALSE, useBytes = FALSE, invert = FALSE){
	x[grep(pattern, x, ignore.case = ignore.case, perl = perl, value = value,
				 fixed = fixed, useBytes = useBytes, invert = invert)]
}

subSampleColumnPerGroup<-function(mat,groupVector,n=NULL){
	lvlTable<-table(groupVector)
	lvls<-names(lvlTable)
	if(is.null(n)){
		n<-min(lvlTable)
	} else {
		if(n > min(lvlTable)) stop("n superior to minimum level freq")
	}
	subSampledData<-lapply(lvls,function(lvl){
		samples <-which(groupVector==lvl)
		selected<-sample(samples,size = n)
		mat[,selected]
	});names(subSampledData)<-lvls
	do.call("cbind",subSampledData)
}

#<Alias
rn<-rownames;
cn<-colnames;
len<-length;
inter<-intersect;
read.tsv<-lire
write.tsv<-ecrire
#>
