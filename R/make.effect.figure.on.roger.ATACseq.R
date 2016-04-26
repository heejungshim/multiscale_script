## `make.effect.figure.on.roger.ATACseq.R' makes effect size figures from multiseq, DESeq for selected sites
#
## Example Usage in /home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/selectedRegions/com R CMD BATCH --no-save --no-restore "--args utils.path='/depot/hjshim/data/hjshim/projects/utils_shim/' hdf5.data.path='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/' library.read.depth.path='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/' sample.prob.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/' pcr.posi.null.path=NULL pcr.posi.alt.path=NULL footprints.path='/home/hjshim/d/shared_data/internal_restricted/roger_atacseq2/dnasefootprints/OpenChromDnaseGm19239.6c.bed.gz' TSS.path='/home/hjshim/d/shared_data/internal_restricted/roger_atacseq2/dnasefootprints/Copper.TSS.DiffExpressed.FDR10.bed.gz' info.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/res/Copper.1024.both.null.f0.all.info' DESeq.100.info.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.null.output.f0/res.Robj' DESeq.300.info.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.300.null.output.f0/res.Robj' DESeq.1024.info.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.1024.null.output.f0/res.Robj' out.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/selectedRegions/' index.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/selectedRegions/Copper.1024.both.null.f0.onlyms.msVSmsmean' file.name='Copper.1024.both.null.f0.onlyms.msVSmsmean' siteSize=1024 treatment='Copper' strand='both' null.permutation=FALSE sig.level=2 multiseq.effect=TRUE deseq.100.effect=TRUE deseq.300.effect=TRUE" /home/hjshim/d/hjshim/projects/multiscale/multiscale_script/R/make.effect.figure.on.roger.ATACseq.R
##
##
## utils.path : path to utils_shim repo ( = '/depot/hjshim/data/hjshim/projects/utils_shim/')
## hdf5.data.path : path to hdf5 file (='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/')
## library.read.depth.path : path to library read depth (='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/')
## sample.prob.path : path to sample prob information (='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/')
## pcr.posi.null.path (default NULL): path to already computed pcr.posi.path for null data 
## pcr.posi.alt.path (default NULL): path to already computed pcr.posi.path for alt data 
## siteSize : site size
## treatment : treatment name
## null.permutation : obtain null using permutation?
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## footprints.path : path to footprints information
## TSS.path : path to TSS information
## info.path : path to file that contains information on regions which are included in our analysis ("index", "chr", "sites", "st.posi", "en.posi", "qval.wave", "qval.ms", "qval.ms.mean", "qval.ms.shape", "qval.deseq.100", "qval.deseq.300", "qval.deseq.1024", "pval.ms", "pval.ms.mean", "pval.ms.shape", "pval.wave", "pval.deseq.100", "pval.deseq.300", "pval.deseq.1024", "logLR.ms.alt", "logLR.ms.null", "logLR.ms.mean.alt", "logLR.ms.mean.null", "logLR.ms.shape.alt", "logLR.ms.shape.null")
## index.path : path to file that contains indexes in info.path file
## DESeq.100.info.path : path to file that contains DESeq2 (bin : 100)  results for sub windows (p-value and fold change). It contains null and alternative data sets together!
## DESeq.300.info.path : path to file that contains DESeq2 (bin : 300)  results for sub windows (p-value and fold change). It contains null and alternative data sets together!
## DESeq.1024.info.path : path to file that contains DESeq2 (bin : 1024) results for sub windows (p-value or fold change). It contains null and alternative data sets together!
## out.path : path to directory where figures will be saved
## file.name : output figure file name
## sig.level : +/- sig.level * standard deviation
## multiseq.effect : indicate whether effect size from multiseq is plotted
## deseq.100.effect : indicate whether effect size from DESeq2 (bin :100) is plotted
## deseq.300.effect : indicate whether effect size from DESeq2 (bin :300) is plotted
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


library("multiseq")
library("ashr")


##utils.path = '/depot/hjshim/data/hjshim/projects/utils_shim/'
##hdf5.data.path ='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/'
##library.read.depth.path='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/'
##sample.prob.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/'
##pcr.posi.null.path=NULL
##pcr.posi.alt.path=NULL
##siteSize=1024
##treatment='Copper'
##strand='both'
##null.permutation=FALSE
##footprints.path='/home/hjshim/d/shared_data/internal_restricted/roger_atacseq2/dnasefootprints/OpenChromDnaseGm19239.6c.bed.gz'
##TSS.path='/home/hjshim/d/shared_data/internal_restricted/roger_atacseq2/dnasefootprints/Copper.TSS.DiffExpressed.FDR10.bed.gz'
##info.path = '/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/res/Copper.1024.both.null.f0.all.info'
##index.path = '/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/selectedRegions/Copper.1024.both.null.f0.both.msVSmsmean'
##DESeq.100.info.path = '/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.100.null.output.f0/res.Robj'
##DESeq.300.info.path = '/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.300.null.output.f0/res.Robj'
##DESeq.1024.info.path = '/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/deseq/Copper.1024.both.1024.null.output.f0/res.Robj'
##out.path = '/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/selectedRegions/'
##file.name= 'Copper.1024.both.null.f0.both.msVSmsmean'
##sig.level = 2
##multiseq.effect=TRUE
##deseq.100.effect=TRUE
##deseq.300.effect=TRUE


args = (commandArgs(TRUE))
eval(parse(text=args))

source(paste0(utils.path, 'R/pcr.artifact.R'))
source(paste0(utils.path, 'R/rhdf5.R'))



## assigen treatment and control name according to input
## treatment    alt     null    control 
## Copper       N702    N705    N706
## Selenium     N703    N705    N706
## Retinoic     N704    N706    N705



##############################################
## sample name and sample file for null data
##############################################

if(null.permutation){
  null = FALSE
}else{
  null = TRUE
}  
name.treatment = NULL
name.control = NULL
if(treatment=='Copper'){
    name.control = "N706"
    if(!null){
        name.treatment = "N702"
    }else{
        name.treatment = "N705"
    }
}
if(treatment=='Selenium'){
    name.control = "N706"
    if(!null){
        name.treatment = "N703"
    }else{
        name.treatment = "N705"
    }
}
if(treatment=='Retinoic'){
    name.control = "N705"
    if(!null){
        name.treatment = "N704"
    }else{
        name.treatment = "N706"
    }
}

## sample names
names.Sam = c("N501", "N502", "N503")

## Make a list of sample names and a list of hdf5 file names : treatment first and control later.
sample.names = c(paste0(name.treatment, names.Sam), paste0(name.control, names.Sam))
sample.files = paste0(sample.names, ".qfiltered10")

sample.names.null = sample.names
sample.files.null = sample.files

##############################################
## sample name and sample file for alternative data
##############################################
null = FALSE
name.treatment = NULL
name.control = NULL
if(treatment=='Copper'){
    name.control = "N706"
    if(!null){
        name.treatment = "N702"
    }else{
        name.treatment = "N705"
    }
}
if(treatment=='Selenium'){
    name.control = "N706"
    if(!null){
        name.treatment = "N703"
    }else{
        name.treatment = "N705"
    }
}
if(treatment=='Retinoic'){
    name.control = "N705"
    if(!null){
        name.treatment = "N704"
    }else{
        name.treatment = "N706"
    }
}

## sample names
names.Sam = c("N501", "N502", "N503")

## Make a list of sample names and a list of hdf5 file names : treatment first and control later.
sample.names = c(paste0(name.treatment, names.Sam), paste0(name.control, names.Sam))
sample.files = paste0(sample.names, ".qfiltered10")

sample.names.alt = sample.names
sample.files.alt = sample.files

numSam = length(sample.names)

## Read sampling prob
if(!is.null(sample.prob.path)){
  if(null.permutation){
    sample.prob.null = scan(paste0(sample.prob.path, treatment, ".alt.prob"), what=double())
  }else{
    sample.prob.null = scan(paste0(sample.prob.path, treatment, ".null.prob"), what=double())
  }
  sample.prob.alt = scan(paste0(sample.prob.path, treatment, ".alt.prob"), what=double())
}else{
  sample.prob.null = sample.prob.alt = rep(1, length(sample.names))
}


## Make a covariate
g = c(rep(0, length(names.Sam)), rep(1, length(names.Sam)))

## pcr posi is given
if(!is.null(pcr.posi.null.path)){
  pcr.posi.null.given = TRUE
}else{
  pcr.posi.null.given = FALSE
}

## pcr posi is given
if(!is.null(pcr.posi.alt.path)){
  pcr.posi.alt.given = TRUE
}else{
  pcr.posi.alt.given = FALSE
}


## read library read depth null
library.read.depth.fwd = rep(0, numSam)
path.read.depth = paste0(library.read.depth.path, "library.read.depth.fwd")
library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
for(i in 1:numSam){
  library.read.depth.fwd[i] = library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names.null[i]),2]
}

library.read.depth.rev = rep(0, numSam)
path.read.depth = paste0(library.read.depth.path, "library.read.depth.rev")
library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
for(i in 1:numSam){
  library.read.depth.rev[i] = library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names.null[i]),2]
}

library.read.depth.fwd.null = library.read.depth.fwd
library.read.depth.rev.null = library.read.depth.rev

## read library read depth alt
library.read.depth.fwd = rep(0, numSam)
path.read.depth = paste0(library.read.depth.path, "library.read.depth.fwd")
library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
for(i in 1:numSam){
  library.read.depth.fwd[i] = library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names.alt[i]),2]
}

library.read.depth.rev = rep(0, numSam)
path.read.depth = paste0(library.read.depth.path, "library.read.depth.rev")
library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
for(i in 1:numSam){
  library.read.depth.rev[i] = library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names.alt[i]),2]
}

library.read.depth.fwd.alt = library.read.depth.fwd
library.read.depth.rev.alt = library.read.depth.rev

## read TF and TSS
all_bed = read.table(gzfile(footprints.path))
TssAnno = read.table(gzfile(TSS.path))


## set up working directory and open figure file
setwd(out.path)
numfig = multiseq.effect + deseq.100.effect + deseq.300.effect + 2
pdf(paste0(out.path, file.name, ".effect.pdf"), width=10, height=7)
nf <- layout(matrix(1:numfig,numfig,1,byrow=TRUE),TRUE)


#############################
# read all information
#############################

index.info = scan(index.path, what=double())
dat.info = read.table(file=info.path, header = TRUE, as.is=TRUE)[index.info,]

if(deseq.100.effect){
  load(DESeq.100.info.path)
  deseq.100.info.pval = deseq.100.info.log2FC = rep(NA, length(use))

  deseq.100.info.pval[use==TRUE] = res$pvalue
  deseq.100.info.log2FC[use==TRUE] = res$log2FoldChange

  window.size = 100  
  numC = siteSize%/%window.size
  
  mat = matrix(deseq.100.info.pval,ncol=numC,byrow=T)
  alt.ix = seq((dim(mat)[1]/2 +1), (dim(mat)[1]), by=1)
  mat.alt = mat[alt.ix,]

  deseq.100.info.pval = mat.alt[dat.info$index,]

  mat = matrix(deseq.100.info.log2FC,ncol=numC,byrow=T)
  alt.ix = seq((dim(mat)[1]/2 +1), (dim(mat)[1]), by=1)
  mat.alt = mat[alt.ix,]

  deseq.100.info.log2FC = mat.alt[dat.info$index,]
}

if(deseq.300.effect){
  load(DESeq.300.info.path)
  deseq.300.info.pval = deseq.300.info.log2FC = rep(NA, length(use))

  deseq.300.info.pval[use==TRUE] = res$pvalue
  deseq.300.info.log2FC[use==TRUE] = res$log2FoldChange

  window.size = 300  
  numC = siteSize%/%window.size
  
  mat = matrix(deseq.300.info.pval,ncol=numC,byrow=T)
  alt.ix = seq((dim(mat)[1]/2 +1), (dim(mat)[1]), by=1)
  mat.alt = mat[alt.ix,]

  deseq.300.info.pval = mat.alt[dat.info$index,]

  mat = matrix(deseq.300.info.log2FC,ncol=numC,byrow=T)
  alt.ix = seq((dim(mat)[1]/2 +1), (dim(mat)[1]), by=1)
  mat.alt = mat[alt.ix,]

  deseq.300.info.log2FC = mat.alt[dat.info$index,]
}

load(DESeq.1024.info.path)
deseq.1024.info.pval = deseq.1024.info.log2FC = rep(NA, length(use))
deseq.1024.info.pval[use==TRUE] = res$pvalue
deseq.1024.info.log2FC[use==TRUE] = res$log2FoldChange

mat = deseq.1024.info.pval
alt.ix = seq((length(mat)/2 +1), (length(mat)), by=1)
mat.alt = mat[alt.ix]
deseq.1024.info.pval = mat.alt[dat.info$index]

mat = deseq.1024.info.log2FC
alt.ix = seq((length(mat)/2 +1), (length(mat)), by=1)
mat.alt = mat[alt.ix]
deseq.1024.info.log2FC = mat.alt[dat.info$index]


numSites = dim(dat.info)[1]

for(ss in 1:numSites){

## ss = 1


sites = dat.info$sites[ss]  
set.seed(sites)
  
## read location information
chr = dat.info$chr[ss]
st.posi = dat.info$st.posi[ss]
en.posi = dat.info$en.posi[ss]

#####################
# read data nulll
#####################
sample.names = sample.names.null
sample.files = sample.files.null
library.read.depth.fwd = library.read.depth.fwd.null
library.read.depth.rev = library.read.depth.rev.null
sample.prob = sample.prob.null
pcr.posi.path = pcr.posi.null.path
pcr.posi.given = pcr.posi.null.given
numSam = length(sample.names)
numBPs = siteSize
library.read.depth = rep(0, numSam)
ATAC.dat = matrix(data=0, nr = numSam, nc = numBPs)
pcr.posi = vector("list", 2)
if(!is.null(pcr.posi.path)){
  pcr.file.path = paste0(pcr.posi.path, "pcrposi.", chr, ".", sites, ".out")
  if(file.exists(pcr.file.path)){
    pcr.posi[[1]] = as.numeric(read.table(file = pcr.file.path, nrows = 1))[-1]
    pcr.posi[[2]] = as.numeric(read.table(file = pcr.file.path, nrows = 1, skip = 1))[-1]
  }else{
    pcr.posi.given = FALSE
  }
}

pcr.ix = 1

## for fwd
if((strand=='both') | (strand=='plus')){

    library.read.depth = library.read.depth + library.read.depth.fwd 

    ## read read counts for a given region
    ## for + strand, we need get reads at locations that are shifted 4bp to left
    ATAC.dat.fwd = matrix(data=NA, nr = numSam, nc = numBPs)
    for(i in 1:numSam){
        path.fwd = paste0(hdf5.data.path, sample.files[i] , ".fwd.h5")
        ATAC.dat.fwd[i, 1:numBPs] = as.matrix(get.counts.h5(path.fwd, paste0("chr", chr), st.posi-4, en.posi-4))
        ## subsample reads..
        if(sample.prob[i] < 1){
          ATAC.dat.fwd[i, 1:numBPs] = rbinom(numBPs, ATAC.dat.fwd[i, 1:numBPs] , sample.prob[i])
        }
    }
    
    ## remove pcr artifacts
    if(!pcr.posi.given){
      pcr.removed.fwd = remove.pcr.artifacts(data=ATAC.dat.fwd, win.half.size=50, prop.thresh=0.9)
      ATAC.dat = ATAC.dat + pcr.removed.fwd$data
      if(!is.null(pcr.removed.fwd$posi.with.pcr.artifacts)){
        pcr.posi[[pcr.ix]] = pcr.removed.fwd$posi.with.pcr.artifacts
      }
    }else{
      pcr.removed.fwd = remove.pcr.artifacts.in.known.posi(data=ATAC.dat.fwd, known.pcr.posi = pcr.posi[[pcr.ix]] , win.half.size = 50, prop.thresh = 0.9)
      ATAC.dat = ATAC.dat + pcr.removed.fwd$data
    }
    pcr.ix = pcr.ix + 1

}

## for reverse
if((strand=='both') | (strand=='minus')){

    library.read.depth = library.read.depth + library.read.depth.rev
  
    ## read read counts for a given region
    ## for - strand, we need get reads at locations that are shifted 4bp to right
    ATAC.dat.rev = matrix(data=NA, nr = numSam, nc = numBPs)
    for(i in 1:numSam){
        path.rev = paste0(hdf5.data.path, sample.files[i] , ".rev.h5")
        ATAC.dat.rev[i, 1:numBPs] = as.matrix(get.counts.h5(path.rev, paste0("chr", chr), st.posi+4, en.posi+4))
        ## subsample reads..
        if(sample.prob[i] < 1){
          ATAC.dat.rev[i, 1:numBPs] = rbinom(numBPs, ATAC.dat.rev[i, 1:numBPs] , sample.prob[i])
        }
    }

    ## remove pcr artifacts
    if(!pcr.posi.given){
      pcr.removed.rev = remove.pcr.artifacts(data=ATAC.dat.rev, win.half.size=50, prop.thresh=0.9)
      ATAC.dat = ATAC.dat + pcr.removed.rev$data
      if(!is.null(pcr.removed.rev$posi.with.pcr.artifacts)){
        pcr.posi[[pcr.ix]] = pcr.removed.rev$posi.with.pcr.artifacts
      }
    }else{
      pcr.removed.rev = remove.pcr.artifacts.in.known.posi(data=ATAC.dat.rev, known.pcr.posi = pcr.posi[[pcr.ix]] , win.half.size=50, prop.thresh=0.9)
      ATAC.dat = ATAC.dat + pcr.removed.rev$data
    } 
    pcr.ix = pcr.ix + 1
    
}

phenoD = ATAC.dat

library.read.depth = library.read.depth*sample.prob 

if(null.permutation){
  ix.permutation = sample(1:6)
  phenoD = phenoD[ix.permutation,]
  library.read.depth = library.read.depth[ix.permutation]
}


####################
# plot raw data
####################
 
## get phenotype
xmin = st.posi
xmax = en.posi

phe.D = phenoD/library.read.depth
trt.pheno = apply(phe.D[1:(numSam/2),], 2, mean)
ctr.pheno = apply(phe.D[(numSam/2+1):numSam,], 2, mean)
trt.RC = sum(phenoD[1:(numSam/2),])
ctr.RC = sum(phenoD[(numSam/2+1):numSam,])


ymin = 0
ymaxT = max(trt.pheno, ctr.pheno)*(1+ 0.05)

xval = xmin:xmax
ymax = ymaxT*10^6


## get pcr information
xval_mapp = NULL
if(length(unlist(pcr.posi)) > 0){ 
    xval_mapp = xval[unlist(pcr.posi)]
}


## Make a raw phenotype figure
raw.title = paste0("chr", chr, ":", st.posi, "-", en.posi, ", control(red):", trt.RC, " control(blue):", ctr.RC)
par(mar = c(1,4,2,2))
plot(1,1,type="n", xlab = "position", ylab = "ATAC-seq rate per million reads",ylim=c(ymin, ymax),xlim=c(xmin, xmax),main =raw.title, axes=FALSE)
##axis(1)
if(!is.null(xval_mapp)){
    axis(1, at=xval_mapp, labels = FALSE, lwd.ticks = 2, col="dark green")
}
axis(2)
box()

### Transcription factor
sel.sites = all_bed[all_bed[,1] == paste("chr", chr, sep="") & all_bed[,2] < (xmax+1) & all_bed[,3] > xmin, ]
offset = -0.0025
if(dim(sel.sites)[1] > 0){
for(k in 1:dim(sel.sites)[1]){
offset = -offset
text(x=(sel.sites[k,2] + sel.sites[k,3])/2, y=(ymax -abs(offset) - offset), strsplit(as.character(sel.sites[k,4]), split="=")[[1]][2])
rect(sel.sites[k,2], 0, sel.sites[k,3], ymax + 1, col=rgb(0,1,0,0.3), border='NA')
}
}


points(xval, ctr.pheno*10^6, col = rgb(0,0,1,alpha=0.7), type="l")
points(xval, trt.pheno*10^6, col = rgb(1,0,0,alpha=0.7), type="l")

#GETS AND PLOTS ANY TSSs IN THE REGION
TSS <- TssAnno[(as.character(TssAnno[,1]) == paste("chr", chr, sep="")) & (TssAnno[,2] > xmin) & (TssAnno[,2] < (xmax+1)),]
if(dim(TSS)[1] > 0) {
for(k in 1:dim(TSS)[1]){
mtext('*', side=1, at=TSS[k,2], col='purple', cex=1.5, padj=1)
}
}


#####################
# read data alt
#####################

sites = dat.info$sites[ss]  
set.seed(sites)

sample.names = sample.names.alt
sample.files = sample.files.alt
library.read.depth.fwd = library.read.depth.fwd.alt
library.read.depth.rev = library.read.depth.rev.alt
sample.prob = sample.prob.alt
pcr.posi.path = pcr.posi.alt.path
pcr.posi.given = pcr.posi.alt.given

numSam = length(sample.names)
numBPs = siteSize
library.read.depth = rep(0, numSam)
ATAC.dat = matrix(data=0, nr = numSam, nc = numBPs)
pcr.posi = vector("list", 2)
if(!is.null(pcr.posi.path)){
  pcr.file.path = paste0(pcr.posi.path, "pcrposi.", chr, ".", sites, ".out")
  if(file.exists(pcr.file.path)){
    pcr.posi[[1]] = as.numeric(read.table(file = pcr.file.path, nrows = 1))[-1]
    pcr.posi[[2]] = as.numeric(read.table(file = pcr.file.path, nrows = 1, skip = 1))[-1]
  }else{
    pcr.posi.given = FALSE
  }
}

pcr.ix = 1

## for fwd
if((strand=='both') | (strand=='plus')){

    library.read.depth = library.read.depth + library.read.depth.fwd 

    ## read read counts for a given region
    ## for + strand, we need get reads at locations that are shifted 4bp to left
    ATAC.dat.fwd = matrix(data=NA, nr = numSam, nc = numBPs)
    for(i in 1:numSam){
        path.fwd = paste0(hdf5.data.path, sample.files[i] , ".fwd.h5")
        ATAC.dat.fwd[i, 1:numBPs] = as.matrix(get.counts.h5(path.fwd, paste0("chr", chr), st.posi-4, en.posi-4))
        ## subsample reads..
        if(sample.prob[i] < 1){
          ATAC.dat.fwd[i, 1:numBPs] = rbinom(numBPs, ATAC.dat.fwd[i, 1:numBPs] , sample.prob[i])
        }
    }
    
    ## remove pcr artifacts
    if(!pcr.posi.given){
      pcr.removed.fwd = remove.pcr.artifacts(data=ATAC.dat.fwd, win.half.size=50, prop.thresh=0.9)
      ATAC.dat = ATAC.dat + pcr.removed.fwd$data
      if(!is.null(pcr.removed.fwd$posi.with.pcr.artifacts)){
        pcr.posi[[pcr.ix]] = pcr.removed.fwd$posi.with.pcr.artifacts
      }
    }else{
      pcr.removed.fwd = remove.pcr.artifacts.in.known.posi(data=ATAC.dat.fwd, known.pcr.posi = pcr.posi[[pcr.ix]] , win.half.size = 50, prop.thresh = 0.9)
      ATAC.dat = ATAC.dat + pcr.removed.fwd$data
    }
    pcr.ix = pcr.ix + 1

}

## for reverse
if((strand=='both') | (strand=='minus')){

    library.read.depth = library.read.depth + library.read.depth.rev
  
    ## read read counts for a given region
    ## for - strand, we need get reads at locations that are shifted 4bp to right
    ATAC.dat.rev = matrix(data=NA, nr = numSam, nc = numBPs)
    for(i in 1:numSam){
        path.rev = paste0(hdf5.data.path, sample.files[i] , ".rev.h5")
        ATAC.dat.rev[i, 1:numBPs] = as.matrix(get.counts.h5(path.rev, paste0("chr", chr), st.posi+4, en.posi+4))
        ## subsample reads..
        if(sample.prob[i] < 1){
          ATAC.dat.rev[i, 1:numBPs] = rbinom(numBPs, ATAC.dat.rev[i, 1:numBPs] , sample.prob[i])
        }
    }

    ## remove pcr artifacts
    if(!pcr.posi.given){
      pcr.removed.rev = remove.pcr.artifacts(data=ATAC.dat.rev, win.half.size=50, prop.thresh=0.9)
      ATAC.dat = ATAC.dat + pcr.removed.rev$data
      if(!is.null(pcr.removed.rev$posi.with.pcr.artifacts)){
        pcr.posi[[pcr.ix]] = pcr.removed.rev$posi.with.pcr.artifacts
      }
    }else{
      pcr.removed.rev = remove.pcr.artifacts.in.known.posi(data=ATAC.dat.rev, known.pcr.posi = pcr.posi[[pcr.ix]] , win.half.size=50, prop.thresh=0.9)
      ATAC.dat = ATAC.dat + pcr.removed.rev$data
    } 
    pcr.ix = pcr.ix + 1
    
}

phenoD = ATAC.dat

library.read.depth = library.read.depth*sample.prob 


####################
# plot raw data
####################

## get phenotype
xmin = st.posi
xmax = en.posi

phe.D = phenoD/library.read.depth
trt.pheno = apply(phe.D[1:(numSam/2),], 2, mean)
ctr.pheno = apply(phe.D[(numSam/2+1):numSam,], 2, mean)
trt.RC = sum(phenoD[1:(numSam/2),])
ctr.RC = sum(phenoD[(numSam/2+1):numSam,])

ymin = 0
ymaxT = max(trt.pheno, ctr.pheno)*(1+ 0.05)

xval = xmin:xmax
ymax = ymaxT*10^6


## get pcr information
xval_mapp = NULL
if(length(unlist(pcr.posi)) > 0){ 
    xval_mapp = xval[unlist(pcr.posi)]
}


## Make a raw phenotype figure
raw.title = paste0("chr", chr, ":", st.posi, "-", en.posi, ", treatment(red):", trt.RC, " control(blue):", ctr.RC)
par(mar = c(1,4,2,2))
plot(1,1,type="n", xlab = "position", ylab = "ATAC-seq rate per million reads",ylim=c(ymin, ymax),xlim=c(xmin, xmax),main =raw.title, axes=FALSE)
axis(1)
if(!is.null(xval_mapp)){
    axis(1, at=xval_mapp, labels = FALSE, lwd.ticks = 2, col="dark green")
}
axis(2)
box()


### Transcription factor
sel.sites = all_bed[all_bed[,1] == paste("chr", chr, sep="") & all_bed[,2] < (xmax+1) & all_bed[,3] > xmin, ]
offset = -0.0025
if(dim(sel.sites)[1] > 0){
for(k in 1:dim(sel.sites)[1]){
offset = -offset
text(x=(sel.sites[k,2] + sel.sites[k,3])/2, y=(ymax -abs(offset) - offset), strsplit(as.character(sel.sites[k,4]), split="=")[[1]][2])
rect(sel.sites[k,2], 0, sel.sites[k,3], ymax + 1, col=rgb(0,1,0,0.3), border='NA')
}
}


points(xval, ctr.pheno*10^6, col = rgb(0,0,1,alpha=0.7), type="l")
points(xval, trt.pheno*10^6, col = rgb(1,0,0,alpha=0.7), type="l")


#GETS AND PLOTS ANY TSSs IN THE REGION
TSS <- TssAnno[(as.character(TssAnno[,1]) == paste("chr", chr, sep="")) & (TssAnno[,2] > xmin) & (TssAnno[,2] < (xmax+1)),]
if(dim(TSS)[1] > 0) {
for(k in 1:dim(TSS)[1]){
mtext('*', side=1, at=TSS[k,2], col='purple', cex=1.5, padj=1)
}
}

########################
# multiseq effect size
########################

if(multiseq.effect){

    ## title
    title = paste0("multiseq [+/-", sig.level, "] -log10(pval): ", round(-log(dat.info$pval.ms[ss],10),2), " logLR: ", round(dat.info$logLR.ms.alt[ss],2), " logLR (shape): ", round(dat.info$logLR.ms.shape.alt[ss],2))
    
    ## get effect size
    genoD = g
    res = multiseq(x = phenoD, g = genoD, read.depth = library.read.depth)

    effect.mean = -res$effect.mean
    effect.sd = sqrt(res$effect.var)

    effect.low = effect.mean - sig.level*effect.sd
    effect.high= effect.mean + sig.level*effect.sd

    ymax = max(effect.high) + 10^(-7)
    ymin = min(effect.low) - 10^(-7)
    
    wh.low = which(effect.low > 0)
    wh.high = which(effect.high < 0)
    high_wh = sort(unique(union(wh.low, wh.high)))
    col_posi = xval[high_wh]

    par(mar = c(1,4,4,2))
    plot(1,1,type="n", xlab = "position", ylab = "Effect size",ylim=c(ymin, ymax),xlim=c(xmin, xmax),main =title, axes=FALSE)
    axis(2)
    if(length(col_posi) > 0){
        for(j in 1:length(col_posi)){
            polygon(c(col_posi[j]-0.5, col_posi[j]-0.5, col_posi[j]+0.5, col_posi[j]+0.5), c(ymin-2, ymax+2, ymax+2, ymin-2), col ="pink", border = NA)
        }
    }
    abline(h = 0, col = "red")
    points(xval, effect.mean, col = "blue", type="l")
    points(xval, effect.high, col = "skyblue", type="l")
    points(xval, effect.low, col = "skyblue", type="l")
    box()

}

########################
## deseq 100 effect
########################

if(deseq.100.effect){

  window.size = 100  
  numC = siteSize%/%window.size

  ## read p-value
  deseq.mlogpval.all = -log(as.numeric(deseq.100.info.pval[ss,]),10)
  deseq.mlogpval.max = max(deseq.mlogpval.all, na.rm = TRUE)
  ## read log2FC
  deseq.log2FC = as.numeric(deseq.100.info.log2FC[ss,])
  
  ## title
  title = paste0("DESeq2 -log10(pval): ", round(-log(dat.info$pval.deseq.100[ss],10),2), " -log(min(pval)): ", round(deseq.mlogpval.max,2), " -log10(pval)[1024]: ", round(-log(dat.info$pval.deseq.1024[ss],10),2))

  ymax.t = 4
  ymin.t = 0

  plot(1,1,type="n", xlab = "position", ylab = "DESeq -log10(pvalue)",ylim=c(ymin.t, ymax.t),xlim=c(xmin, xmax),main = title, axes=FALSE)
  axis(2)
  xleft = rep(NA, numC)
  xright = rep(NA, numC)
  xleft[1] = xmin
  for(j in 1:(numC-1)){
    xleft[j+1] = xleft[j] + window.size
    xright[j] = xleft[j+1] - 1
  }
  xright[numC] = xmax
    
  ybottom = ytop = rep(0,numC)
  ytop = deseq.mlogpval.all
  wh = which(is.na(ytop)==TRUE)
  if(length(wh) > 0){
    ytop[wh] = 0
  }
  rect(xleft, ybottom, xright, ytop, col = "grey")
  box()

}




########################
## deseq 300 effect
########################

if(deseq.300.effect){

  window.size = 300  
  numC = siteSize%/%window.size

  ## read p-value
  deseq.mlogpval.all = -log(as.numeric(deseq.300.info.pval[ss,]),10)
  deseq.mlogpval.max = max(deseq.mlogpval.all, na.rm = TRUE)
  ## read log2FC
  deseq.log2FC = as.numeric(deseq.300.info.log2FC[ss,])

  ## title
  title = paste0("DESeq2 -log10(pval): ", round(-log(dat.info$pval.deseq.300[ss],10),2), " -log(min(pval)): ", round(deseq.mlogpval.max,2), " logLR (mean): ", round(dat.info$logLR.ms.mean.alt[ss],2))

  ymax.t = 4
  ymin.t = 0

  plot(1,1,type="n", xlab = "position", ylab = "DESeq -log10(pvalue)",ylim=c(ymin.t, ymax.t),xlim=c(xmin, xmax),main = title, axes=FALSE)
  axis(2)
  xleft = rep(NA, numC)
  xright = rep(NA, numC)
  xleft[1] = xmin
  for(j in 1:(numC-1)){
    xleft[j+1] = xleft[j] + window.size
    xright[j] = xleft[j+1] - 1
  }
  xright[numC] = xmax
    
  ybottom = ytop = rep(0,numC)
  ytop = deseq.mlogpval.all
  wh = which(is.na(ytop)==TRUE)
  if(length(wh) > 0){
    ytop[wh] = 0
  }  
  rect(xleft, ybottom, xright, ytop, col = "grey")
  box()
  
}

}

dev.off()


 
