## `preprocess.wave.DESeq.N.run.multiscale.on.roger.ATACseq.R' prepare data for WaveQTL, DESeq2, and run multiseq 
##
##
## Example Usage (see command in /home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/gen.data/com/Copper.1024.both.alt) :
## R CMD BATCH --no-save --no-restore "--args chr=1 sites.ix=1 sites.iv=1000 pcr.posi.path=NULL pcr.posi.print=TRUE com.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/gen.data/com/' WaveQTL.path='/home/hjshim/d/hjshim/software/WaveQTL/' utils.path='/depot/hjshim/data/hjshim/projects/utils_shim/' hdf5.data.path='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/' library.read.depth.path='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/' loc.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/locus/' sample.prob.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/' wd.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/run/' siteSize=1024 treatment='Copper' strand='both' null=FALSE null.permutation=FALSE meanR.thresh=2 window.size.list=c(100,300,1024) wavelet.preprocess.QT=TRUE wavelet.preprocess.NoQT=FALSE deseq.preprocess=TRUE" /home/hjshim/d/hjshim/projects/multiscale/multiscale_script/R/preprocess.wave.DESeq.N.run.multiscale.on.roger.ATACseq.R
##
##
## WaveQTL.path : path to WaveQTL software (= '/home/hjshim/d/hjshim/software/WaveQTL/')
## utils.path : path to utils_shim repo ( = '/depot/hjshim/data/hjshim/projects/utils_shim/')
## hdf5.data.path : path to hdf5 file (='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/')
## library.read.depth.path : path to library read depth (='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/')
## loc.path : path to locus informatin (='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/locus/')
## sample.prob.path : path to sample prob information (='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/')
## wd.path: working directory path ('/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/run/')
## chr : chromosome name  
## sites.ix : see sites.iv 
## sites.iv : with sites.ix, we will define st.sties and en.sitnes (st.sites = sites.iv * (sites.ix -1) + 1 en.sites = min(sites.iv * sites.ix, numSites)) and we will run from st.sites to en.sites.
## pcr.posi.path (default NULL): path to already computed pcr.posi.path 
## pcr.posi.print (= TRUE): print out pcr.posi information?
## siteSize : site size
## treatment : treatment name
## null : indicate whether it's null (control 1 vs control 2) or alternative data
## null.permutation : obtain null using permutation?
## strand : 'both', 'plus', 'minus'; add two strands, use + strand, or use - strand
## meanR.thresh=1 : for wavelet preprocess
## window.size.list=c(100,300) : list of window sizes we consider for DESeq analysis
## wavelet.preprocess.QT=TRUE
## wavelet.preprocess.NoQT=TRUE
## deseq.preprocess=TRUE
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

##WaveQTL.path = '/home/hjshim/d/hjshim/software/WaveQTL/'
##utils.path = '/depot/hjshim/data/hjshim/projects/utils_shim/'
##hdf5.data.path ='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/'
##library.read.depth.path='/depot/hjshim/data/shared_data/internal_restricted/roger_atacseq2/hdf5/'
##loc.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/locus/'
##sample.prob.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/info/'
##wd.path='/depot/hjshim/data/hjshim/projects/multiscale/atacseq_analysis/run/' 
##chr=22
##sites.ix=2
##sites.iv=100
##pcr.posi.path=NULL
##pcr.posi.print = TRUE

##siteSize=1024
##treatment='Copper'
##strand='both'
##null=TRUE
##null.permutation=FALSE
##meanR.thresh=2
##window.size.list=c(100,300,1024)
##wavelet.preprocess.QT=TRUE
##wavelet.preprocess.NoQT=FALSE
##deseq.preprocess=TRUE

args = (commandArgs(TRUE))
eval(parse(text=args))

source(paste0(utils.path, 'R/pcr.artifact.R'))
source(paste0(utils.path, 'R/rhdf5.R'))

## assigen treatment and control name according to input
## treatment    alt     null    control 
## Copper       N702    N705    N706
## Selenium     N703    N705    N706
## Retinoic     N704    N706    N705

if(null.permutation){
  null=FALSE
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

## Read sampling prob
if(!is.null(sample.prob.path)){
  if(null){
    sample.prob = scan(paste0(sample.prob.path, treatment, ".null.prob"), what=double())
  }else{
    sample.prob = scan(paste0(sample.prob.path, treatment, ".alt.prob"), what=double())
  }
}else{
  sample.prob = rep(1, length(sample.names))
}

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
if(null.permutation){
  output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".nullpermutation")
}else{
  if(!null){
    output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".alt")
  }else{
    output.dir.name=paste0(treatment,".", siteSize, ".", strand, ".null")
  }
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

## make output directory name and output directory for wavelets
if(wavelet.preprocess.QT){  
    wave.out.dir.path = paste0(wd.path, "wave/", output.dir.name, ".data/") 
    if(!file.exists(wave.out.dir.path)){
        dir.create(wave.out.dir.path)
    }
}

if(wavelet.preprocess.NoQT){  
    waveNoQT.out.dir.path = paste0(wd.path, "waveNoQT/", output.dir.name, ".data/") 
    if(!file.exists(waveNoQT.out.dir.path)){
        dir.create(waveNoQT.out.dir.path)
    }
}

## make output directory name and output directory for DESeq
if(deseq.preprocess){
    for(ss in 1:length(window.size.list)){
      window.size = window.size.list[ss]
      if(null.permutation){
        output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".nullpermutation")
      }else{
        if(!null){
          output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".alt")
        }else{
          output.dir.name=paste0(treatment,".", siteSize, ".", strand, ".", window.size, ".null")
        }
      }
      
      deseq.out.dir.path = paste0(wd.path, "deseq/", output.dir.name, ".data/") 
      if(!file.exists(deseq.out.dir.path)){
        dir.create(deseq.out.dir.path)
      }
    }
  }    


#############################
# read location information and numSites for each chromosome
#############################
path = paste0(loc.path, treatment, ".", siteSize, ".chr", chr, ".locus")
loc.info = read.table(path, header=TRUE)

path = paste0(loc.path, treatment, ".", siteSize, ".numSites.txt")
numSites = scan(path)[chr]


st.sites = sites.iv * (sites.ix -1) + 1
en.sites = min(sites.iv * sites.ix, numSites) 

if(st.sites <= numSites){
for(sites in st.sites:en.sites){

set.seed(sites)
  
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
  ## Heejung working on here!
  ix.permutation = sample(1:6)
  phenoD = phenoD[ix.permutation,]
  library.read.depth = library.read.depth[ix.permutation]
}


#########################################
# data preprocessing for wavelet analysis
#########################################

if(wavelet.preprocess.QT){

    source(paste0(WaveQTL.path, "/R/WaveQTL_preprocess_funcs.R"))
    res = WaveQTL_preprocess(Data = phenoD, library.read.depth = library.read.depth , Covariates = NULL, meanR.thresh = meanR.thresh)

    filteredWCs = res$filtered.WCs
    norm.WCs = res$WCs

    this.path = paste0(wave.out.dir.path, "WC.", chr, ".", sites, ".txt")
    write.table(norm.WCs, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    this.path = paste0(wave.out.dir.path, "use.", chr, ".", sites, ".txt")
    cat(filteredWCs, file = this.path)
  }

if(wavelet.preprocess.NoQT){

    source(paste0(WaveQTL.path, "/R/WaveQTL_preprocess_funcs.R"))
    res = WaveQTL_preprocess(Data = phenoD, library.read.depth = library.read.depth , Covariates = NULL, meanR.thresh = meanR.thresh, no.QT = TRUE)

    filteredWCs = res$filtered.WCs
    norm.WCs = res$WCs

    this.path = paste0(waveNoQT.out.dir.path, "WC.", chr, ".", sites, ".txt")
    write.table(norm.WCs, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    this.path = paste0(waveNoQT.out.dir.path, "use.", chr, ".", sites, ".txt")
    cat(filteredWCs, file = this.path)
}

 
#########################################
# data preprocessing for DESeq
#########################################
if(deseq.preprocess){

    for(ss in 1:length(window.size.list)){
        window.size = window.size.list[ss]
        numC = numBPs%/%window.size
        mat = matrix(data=NA, nc = numSam, nr = numC)
        st = 1
        if(numC > 1){
          for(c in 1:(numC-1)){
            en = st + window.size  - 1
            mat[c,] = apply(phenoD[,st:en], 1, sum)
            st = en + 1
          }
        }
        en = numBPs
        mat[numC,] = apply(phenoD[,st:en], 1, sum)
  
        ## save deseq data
        if(null.permutation){
          output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".nullpermutation")
        }else{
          if(!null){
            output.dir.name=paste0(treatment,".", siteSize, ".", strand,  ".", window.size, ".alt")
          }else{
            output.dir.name=paste0(treatment,".", siteSize, ".", strand, ".", window.size, ".null")
          }
        }

        deseq.out.dir.path = paste0(wd.path, "deseq/", output.dir.name, ".data/")
        this.path = paste0(deseq.out.dir.path, "data.", chr, ".", sites, ".txt")
        write.table(mat, file = this.path, quote= FALSE, row.names = FALSE, col.names = FALSE)
    }

}


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

