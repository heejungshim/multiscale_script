## `run.multiscale.on.roger.ATACseq.R' runs only multiseq on roger ATAC-seq data. It is useful when we try to rerun only multiseq after preprocessing (there is an option to avoid pcr.posi computation).
##
##
## Example Usage (in /mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/multiscale/com/Copper.1024.both.null/) R CMD BATCH --no-save --no-restore "--args chr=10 sites.ix=$SGE_TASK_ID wd.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/' siteSize=1024 treatment='Copper' null=TRUE strand='both' pcr.posi.path='/mnt/lustre/home/shim/multiscale_analysis/analysis/roger_ATAC2/run/multiscale.backup/Copper.1024.both.null.output/' pcr.posi.print=TRUE" /mnt/lustre/home/shim/multiscale_analysis/src/R/run.multiscale.on.roger.ATACseq.R
##
##
## chr : chromosome
## sites.ix : default=NULL; if it is null, run multiseq on all sites or we can specifiy partifular site
## wd.path : working directory path
## siteSize : site size
## treatment : treatment name
## null : indicate whether it's null (control 1 vs control 2) or alternative data
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## pcr.posi.path : default = NULL; if it is null, recompute pcr posi. Otherwise, read pre-computed pcr.posi
## pcr.posi.print : default = TRUE; whether to print pcr posi or not.
##
## Copyright (C) 2015 Heejung Shim
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


##pcr.artifact.script.path='/depot/hjshim/data/hjshim/projects/utils_shim/R/pcr.artifact.R'
##rhdf5.script.path='/depot/hjshim/data/hjshim/projects/utils_shim/R/rhdf5.R'
##hdf5.data.path ='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/'
##library.read.depth.path='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/'
##loc.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/locus/'

##wd.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/run/' 
##siteSize=1024
##treatment='Copper'
##strand='both'
##null=FALSE

##chr=1 
##sites.ix=1
##sites.iv=1
##pcr.posi.path=NULL
##pcr.posi.print = TRUE


args = (commandArgs(TRUE))
eval(parse(text=args))
##nargs()

source(pcr.artifact.script.path)
source(rhdf5.script.path)

## assigen treatment and control name according to input
## treatment    alt     null    control 
## Copper       N702    N705    N706
## Selenium     N703    N705    N706
## Retinoic     N704    N706    N705

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

## Make a covariate
g = c(rep(0, length(names.Sam)), rep(1, length(names.Sam)))

## pcr posi is given
if(!is.null(pcr.posi.path)){
  pcr.posi.given = TRUE
}else{
  pcr.posi.given = FALSE
}

## set up working directory 
setwd(wd.path)

## make output directory name and output directory
if(!null){
    output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".alt")
}else{
    output.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
}


multiseq.out.dir.path = paste0(wd.path, "multiscale/", output.dir.name, ".output") 
if(!file.exists(multiseq.out.dir.path)){
    dir.create(multiseq.out.dir.path)
}


## make warning message directory
warning.dir.path = paste0(wd.path, "multiscale/", output.dir.name, ".warnings") 
if(!file.exists(warning.dir.path)){
    dir.create(warning.dir.path)
}

#############################
## read location information and numSites for each chromosome
#############################
path = paste0(loc.path, treatment, ".", siteSize, ".chr", chr, ".locus")
loc.info = read.table(path, header=TRUE)

path = paste0(loc.path, treatment, ".", siteSize, ".numSites.txt")
numSites = scan(path)[chr]


st.sites = sites.iv * (sites.ix -1) + 1
en.sites = min(sites.iv * sites.ix, numSites) 

if(st.sites <= numSites){
  
for(sites in st.sites:en.sites){

## sites = 1    
## create file for warning message
warn.path = paste0(warning.dir.path, "/warnings.", chr, ".", sites, ".txt")
warnings.file <- file(warn.path, open="wt")
sink(warnings.file, type="message")



#############################
# read location information 
#############################
st.posi = loc.info$st.posi[sites]
en.posi = loc.info$en.posi[sites]


#####################
# read data
#####################
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

    ## read library read depth
    path.read.depth = paste0(library.read.depth.path, "library.read.depth.fwd")
    library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
    for(i in 1:numSam){
        library.read.depth[i] = library.read.depth[i] + library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names[i]),2]
    }

    ## read read counts for a given region
    ## for + strand, we need get reads at locations that are shifted 4bp to left
    ATAC.dat.fwd = matrix(data=NA, nr = numSam, nc = numBPs)
    for(i in 1:numSam){
        path.fwd = paste0(hdf5.data.path, sample.files[i] , ".fwd.h5")
        ATAC.dat.fwd[i, 1:numBPs] = as.matrix(get.counts.h5(path.fwd, paste0("chr", chr), st.posi-4, en.posi-4))
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

    ## read library read depth
    path.read.depth = paste0(library.read.depth.path, "library.read.depth.rev")
    library.read.depth.dat = read.table(path.read.depth, as.is=TRUE)
    for(i in 1:6){
        library.read.depth[i] = library.read.depth[i] + library.read.depth.dat[which(library.read.depth.dat[,1] == sample.names[i]),2]
    }
    
    ## read read counts for a given region
    ## for - strand, we need get reads at locations that are shifted 4bp to right
    ATAC.dat.rev = matrix(data=NA, nr = numSam, nc = numBPs)
    for(i in 1:numSam){
        path.rev = paste0(hdf5.data.path, sample.files[i] , ".rev.h5")
        ATAC.dat.rev[i, 1:numBPs] = as.matrix(get.counts.h5(path.rev, paste0("chr", chr), st.posi+4, en.posi+4))
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


## perform test
genoD = g
res = multiseq(x = phenoD, g = genoD, read.depth = library.read.depth, verbose = TRUE)
out.res = c(res$logLR$value, res$logLR$scales)

## write output
write.table(t(out.res), file = paste0(multiseq.out.dir.path, "/res.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)
if(pcr.posi.print){
  for(m in 1:(pcr.ix-1)){
    if(m == 1){
      write.table(t(c(m, pcr.posi[[m]])), file = paste0(multiseq.out.dir.path, "/pcrposi.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE)
    }else{
      write.table(t(c(m, pcr.posi[[m]])), file = paste0(multiseq.out.dir.path, "/pcrposi.", chr, ".", sites, ".out"), quote= FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}

}
}
