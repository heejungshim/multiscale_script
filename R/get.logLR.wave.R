## `get.logLR.wave.R' contains scrits to collect logLR from wavelet analysis, save them as a vector, and output them as R object.
##
## Example Usage (see command in /home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/wave/com/Copper.1024.both.null.post) : R CMD BATCH --no-save --no-restore "--args chr=1 st.sites=1 en.sites=28978 path.output.dir='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/wave/Copper.1024.both.null.run/' output.file.name='Copper.1024.both.null'" /home/hjshim/d/hjshim/projects/multiscale/multiscale_script/R/get.logLR.wave.R
##
##
## chr : chromosome (if chr=NULL, we assume there is no chromosome). Otherwise, output file format is chr.sites
## st.sites : first site to combine
## en.sites :  last site to combine 
## path.output.dir : path to directory which contain output directory
## output.file.name : output file name in output directory
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



##chr=1
##st.sites=1
##en.sites=100
##path.output.dir='/home/hjshim/d/hjshim/projects/multiscale/atacseq_analysis/run/wave/Copper.1024.both.alt.run/'
##output.file.name='Copper.1024.both.alt'


args = (commandArgs(TRUE))
eval(parse(text=args))  

numSites = en.sites - st.sites + 1
done = rep(NA, numSites)
logLR = rep(NA, numSites)

path = paste0(path.output.dir, "output/", output.file.name, ".")

for(i in st.sites:en.sites){

  if(is.null(chr)){
    path.each  = paste0(path, i, ".fph.logLR.txt")
  }else{
    path.each  = paste0(path, chr, ".", i, ".fph.logLR.txt")
  }
  
  if(file.exists(path.each)== FALSE){		
    done[i] = FALSE
  }else{
    if(file.info(path.each)$size == 0){
      done[i] = FALSE
    }else{			
      dat = scan(path.each, what="")
      done[i] = TRUE
      if(!is.na(as.numeric(dat[2]))){
        logLR[i] = as.numeric(dat[2])
      }
    }
  }
}

## make summary directory
path.sum.dir = paste0(path.output.dir, "sum/")
if(!file.exists(path.sum.dir)){
    dir.create(path.sum.dir)
}

if(is.null(chr)){
    save("logLR", "done", file =paste0(path.sum.dir, output.file.name, ".Robj"))
}else{
    save("logLR", "done", file =paste0(path.sum.dir, output.file.name, ".", chr, ".Robj"))
}






