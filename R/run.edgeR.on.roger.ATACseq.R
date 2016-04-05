## `run.edgeR.on.roger.ATACseq.R' contains scrits to run edgeR on Roger's ATAC data.
## 
##
## Example Usage (see command in /home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/edgeR/Copper.1024.both.100.null.output.f0/com) : R CMD BATCH --no-save --no-restore "--args wd.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/edgeR/' alt.data.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.alt.run/' null.data.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.null.run/' output.dir.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/edgeR/Copper.1024.both.100.null.output.f0/' numSites.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/locus/Copper.1024.numSites.txt' library.read.depth.fwd.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/library.read.depth.fwd.Copper.alt' library.read.depth.rev.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/library.read.depth.rev.Copper.alt' siteSize=1024 window.size=100 numSam=6 filter.cut=0" /home/hjshim/d/hjshim/projects/multiscale/multiscale <- script/R/run.edgeR.on.roger.ATACseq.R
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


##wd.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/edgeR/'
##alt.data.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.alt.run/'
##null.data.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.null.run/'
##output.dir.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/edgeR/Copper.1024.both.100.output.f0/'
##numSites.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/locus/Copper.1024.numSites.txt'
##siteSize=1024
##window.size=100
##numSam = 6
##filter.cut = 0

##library.read.depth.fwd.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/library.read.depth.fwd.Copper.alt'
##library.read.depth.rev.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/library.read.depth.rev.Copper.alt'

args = (commandArgs(TRUE))
eval(parse(text=args))

library("edgeR")

## set up working directory 
setwd(wd.path)

## get the library read depth 
library.read.depth.fwd = scan(library.read.depth.fwd.path, what=double())
library.read.depth.rev = scan(library.read.depth.rev.path, what=double())
library.read.depth = library.read.depth.fwd + library.read.depth.rev

## get number of sites informaiton 
numSites.list = scan(numSites.path)

## read data 
input.dir.path = alt.data.path
numC = siteSize%/%window.size
edgeR.data = matrix(data=NA, nc = numSam, nr = sum(numSites.list)*numC)
st.chr.ix = NULL
en.chr.ix = NULL
st.chr.ix[1] = 1
for(chr in 1:22){
  numSites = numSites.list[chr]
  numRow = numSites*numC
  en.chr.ix[chr] = st.chr.ix[chr] + numRow - 1
  edgeR.data.each = read.table(paste0(input.dir.path, "data.", chr, ".txt"))
  if(dim(edgeR.data.each)[1] == numRow){
    edgeR.data[st.chr.ix[chr]:en.chr.ix[chr],] = as.matrix(edgeR.data.each)
    }
  st.chr.ix[chr+1] = en.chr.ix[chr] + 1
}
alt.data = edgeR.data

## read null data if null.data.path is provided 
if(!is.null(null.data.path)){
  input.dir.path = null.data.path
  numC = siteSize%/%window.size
  edgeR.data = matrix(data=NA, nc = numSam, nr = sum(numSites.list)*numC)
  st.chr.ix = NULL
  en.chr.ix = NULL
  st.chr.ix[1] = 1
  for(chr in 1:22){
    numSites = numSites.list[chr]
    numRow = numSites*numC
    en.chr.ix[chr] = st.chr.ix[chr] + numRow - 1
    edgeR.data.each = read.table(paste0(input.dir.path, "data.", chr, ".txt"))
    if(dim(edgeR.data.each)[1] == numRow){
      edgeR.data[st.chr.ix[chr]:en.chr.ix[chr],] = as.matrix(edgeR.data.each)
    }
    st.chr.ix[chr+1] = en.chr.ix[chr] + 1
  }
  null.data = edgeR.data
}

if(!is.null(null.data.path)){
  edgeR.data = rbind(null.data, alt.data)
}else{
  edgeR.data = alt.data
}

## filter data
rsum = rowSums(edgeR.data)
use = ((rsum > filter.cut) & (!is.na(rsum)))
countData.filtered = edgeR.data[ use, ]

## Perform edgeR Analysis
gp = factor(c(rep(0,3), rep(1,3)) ,labels=c("A","B"))
res = DGEList(counts=countData.filtered, group=gp, lib.size=library.read.depth)
res = calcNormFactors(res, method="RLE")
res = estimateCommonDisp(res)
res = estimateTagwiseDisp(res)
res = exactTest(res, dispersion="auto")

## get p-value
pval.vec = rep(NA, length(use))
pval.vec[use==TRUE] = res$table$PValue
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
