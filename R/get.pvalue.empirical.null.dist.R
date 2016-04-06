## `get.pvalue.empirical.null.dist.R' computes pvalue using 'get.pval.from.empirical.null.dist.discrete' in 'utils_shim' repo.
##
## Example Usage (see command in /home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/etc/Copper.1024.both/com) R CMD BATCH --no-save --no-restore "--args wd.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/etc/Copper.1024.both/' utils.path='/depot/hjshim/data/hjshim/projects/utils_shim/' method='ms' seed=1" /home/hjshim/d/hjshim/projects/multiscale/multiscale <- script/R/get.pvalue.empirical.null.dist.R
##
##
## wd.path: path to working directory 
## utils.path: path to utils_shim repo.
## method (e.g., 'ms'): pvalue from which method?
##
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


##wd.path='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/etc/Copper.1024.both/'
##utils.path = '/depot/hjshim/data/hjshim/projects/utils_shim/'
##method='ms'
##seed = 1

args = (commandArgs(TRUE))
eval(parse(text=args))  

set.seed(seed)

load(paste0(wd.path, "test.stat.Robj"))
## [1] "all.name"     "deseq.100.a"  "deseq.100.n"  "deseq.1024.a" "deseq.1024.n"
## [6] "deseq.300.a"  "deseq.300.n"  "dir.path"     "edgeR.100.a"  "edgeR.100.n" 
##[11] "edgeR.1024.a" "edgeR.1024.n" "edgeR.300.a"  "edgeR.300.n"  "in.ix"       
##[16] "ms.a"         "ms.mean.a"    "ms.mean.n"    "ms.n"         "ms.shape.a"  
##[21] "ms.shape.n"   "numT"         "wave.a"       "wave.n"     
source(paste0(utils.path, "R/permutation.R"))

if(method=='wave'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = wave.n, statistic.alt = wave.a)
  out.dir.path = paste0(dir.path, "wave/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".Robj"))
}

if(method=='ms'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = ms.n, statistic.alt = ms.a)
  out.dir.path = paste0(dir.path, "multiscale/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".Robj"))
}
if(method=='ms.shape'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = ms.shape.n, statistic.alt = ms.shape.a)
  out.dir.path = paste0(dir.path, "multiscale/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".shape.Robj"))
}
if(method=='ms.mean'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = ms.mean.n, statistic.alt = ms.mean.a)
  out.dir.path = paste0(dir.path, "multiscale/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".mean.Robj"))
}
if(method=='deseq.100'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = deseq.100.n, statistic.alt = deseq.100.a, big.sig = FALSE)
  out.dir.path = paste0(dir.path, "deseq/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".100.Robj"))
}
if(method=='deseq.300'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = deseq.300.n, statistic.alt = deseq.300.a, big.sig = FALSE)
  out.dir.path = paste0(dir.path, "deseq/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".300.Robj"))
}
if(method=='deseq.1024'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = deseq.1024.n, statistic.alt = deseq.1024.a, big.sig = FALSE)
  out.dir.path = paste0(dir.path, "deseq/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".1024.Robj"))
}
if(method=='edgeR.100'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = edgeR.100.n, statistic.alt = edgeR.100.a, big.sig = FALSE)
  out.dir.path = paste0(dir.path, "edgeR/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".100.Robj"))
}
if(method=='edgeR.300'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = edgeR.300.n, statistic.alt = edgeR.300.a, big.sig = FALSE)
  out.dir.path = paste0(dir.path, "edgeR/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".300.Robj"))
}
if(method=='edgeR.1024'){
  pval.list = get.pval.from.empirical.null.dist.discrete(statistic.null = edgeR.1024.n, statistic.alt = edgeR.1024.a, big.sig = FALSE)
  out.dir.path = paste0(dir.path, "edgeR/summary/")
  if(!file.exists(out.dir.path)){
    dir.create(out.dir.path)
  }
  pval = rep(NA, numT)
  pval[in.ix] = pval.list
  save("pval", file =paste0(out.dir.path, all.name, ".1024.Robj"))
}




