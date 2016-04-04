## `run.DESeq2.on.roger.ATACseq.R' contains scrits to run DESeq2 on Roger's ATAC data.
## 
##
## Example Usage (see command in /home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.null.output.f0/com) : R CMD BATCH --no-save --no-restore "--args wd.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/' alt.data.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.alt.run/' null.data.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.null.run/' output.dir.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.null.output.f0/' numSites.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/locus/Copper.1024..numSites.txt' siteSize=1024 window.size=100 numSam=6 filter.cut=0" /home/hjshim/d/hjshim/projects/multiscale/multiscale_script/R/run.DESeq2.on.roger.ATACseq.R
##
##
## wd.path : working directory path
## siteSize : site size
## treatment : treatment name
## null : indicate whether it's null (control 1 vs control 2) or alternative data
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## window.size : window size we consider for DESeq analysis
## numSam : number of samples
## filter.cut : analysis includes window with read count > filter.cut
##
##
## Copyright (C) 2014 Heejung Shim
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.


##wd.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/'
##alt.data.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.alt.run/'
##null.data.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.null.run/'
##output.dir.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.output.f0/'
##numSites.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/locus/Copper.1024.numSites.txt'
##siteSize=1024
##window.size=100
##numSam = 6
##filter.cut = 0

args = (commandArgs(TRUE))
eval(parse(text=args))

library("DESeq2")

## set up working directory 
setwd(wd.path)

## set up name list to run DESeq
genoD = c(rep("T", numSam/2), rep("C", numSam/2))
name.list = genoD
for(i in 1:floor(numSam/2)){
  name.list[i] = paste0(genoD[i], i)
  name.list[numSam/2 + i] = paste0(genoD[numSam/2 + i], i)
}


## get number of sites informaiton 
numSites.list = scan(numSites.path)

## read data 
input.dir.path = alt.data.path
numC = siteSize%/%window.size
deseq.data = matrix(data=NA, nc = numSam, nr = sum(numSites.list)*numC)
st.chr.ix = NULL
en.chr.ix = NULL
st.chr.ix[1] = 1
for(chr in 1:22){
  numSites = numSites.list[chr]
  numRow = numSites*numC
  en.chr.ix[chr] = st.chr.ix[chr] + numRow - 1
  deseq.data.each = read.table(paste0(input.dir.path, "data.", chr, ".txt"))
  if(dim(deseq.data.each)[1] == numRow){
    deseq.data[st.chr.ix[chr]:en.chr.ix[chr],] = as.matrix(deseq.data.each)
    }
  st.chr.ix[chr+1] = en.chr.ix[chr] + 1
}
alt.data = deseq.data

## read null data if null.data.path is provided 
if(!is.null(null.data.path)){
  input.dir.path = null.data.path
  numC = siteSize%/%window.size
  deseq.data = matrix(data=NA, nc = numSam, nr = sum(numSites.list)*numC)
  st.chr.ix = NULL
  en.chr.ix = NULL
  st.chr.ix[1] = 1
  for(chr in 1:22){
    numSites = numSites.list[chr]
    numRow = numSites*numC
    en.chr.ix[chr] = st.chr.ix[chr] + numRow - 1
    deseq.data.each = read.table(paste0(input.dir.path, "data.", chr, ".txt"))
    if(dim(deseq.data.each)[1] == numRow){
      deseq.data[st.chr.ix[chr]:en.chr.ix[chr],] = as.matrix(deseq.data.each)
    }
    st.chr.ix[chr+1] = en.chr.ix[chr] + 1
  }
  null.data = deseq.data
}

if(!is.null(null.data.path)){
  deseq.data = rbind(null.data, alt.data)
}else{
  deseq.data = alt.data
}


## Converts the colData table to a dataframe with vector labels
colData <- data.frame(row.names = as.vector(name.list), t = as.vector(genoD), r=as.factor(c(1,2,3,1,2,3)))

## filter data
rsum = rowSums ( deseq.data)
use = ((rsum > filter.cut) & (!is.na(rsum)))
countData.filtered = deseq.data[ use, ]

## Perform DESeq2 Analysis
ddsTvC <- DESeqDataSetFromMatrix(
	countData=countData.filtered,
	colData=colData,
	design = ~ r + t)
dds <- DESeq(ddsTvC)
res <- results(dds)

## get p-value
pval.vec = rep(NA, length(use))
pval.vec[use==TRUE] = res$pvalue
pval.filtered =matrix(pval.vec,ncol=numC,byrow=T)

## get minimum p-value for each site
min.pval=apply(pval.filtered,1,min,na.rm=TRUE)
min.pval[is.infinite(min.pval)] = NA

## try to save output
## make an output directory
if(!file.exists(output.dir.path)){
  dir.create(output.dir.path)
}

## output min.pval
write.table(min.pval, file = paste0(output.dir.path, "/min.pval.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)

## output object just in case we want to take a close look!
save("res", "use", file =paste0(output.dir.path, "/res.Robj"))

## When window size is the same as siteSize, we will compute logLL, logFC, selogFC. 
if(window.size == siteSize){

  mlogLR = qchisq(1-res@listData$pvalue, 1)/2
  logFC = log(2)*res@listData$log2FoldChange
  selogFC = log(2)*res@listData$lfcSE

  mlogLR.vec = logFC.vec = selogFC.vec = rep(NA, length(use))
  mlogLR.vec[use==TRUE] = mlogLR
  logFC.vec[use==TRUE] = logFC
  selogFC.vec[use==TRUE] = selogFC

  write.table(mlogLR.vec, file = paste0(output.dir.path, "/logLR.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
  write.table(logFC.vec, file = paste0(output.dir.path, "/logFC.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
  write.table(selogFC.vec, file = paste0(output.dir.path, "/selogFC.txt"), quote= FALSE, row.names = FALSE, col.names = FALSE)
}












