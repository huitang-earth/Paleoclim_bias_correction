# how to read negative values (signed data) from the raster files?
# how to use R to calculate bioclimate variables

library(rgdal)
library(raster)
library(dismo)
library(ncdf4)

# source("/media/LaCie/linux_backup/home/Hui/cclm/Postprocess/rscripts/fun.gaussian.filter.r")

jobname<-c("11001tort_echam5")

setwd(paste0("/home/huitang/Documents/data/Indre/Esther_redescription_paleo/",jobname))
  
fprec<-nc_open("wc2.0_10m_prec_plio.nc")
prec<-ncvar_get(fprec,"Band1")
ftavg<-nc_open("wc2.0_10m_tavg_plio.nc")
tavg<-ncvar_get(ftavg,"Band1")
ftmax<-nc_open("wc2.0_10m_tmax_plio.nc")
tmax<-ncvar_get(ftmax,"Band1")
ftmin<-nc_open("wc2.0_10m_tmin_plio.nc")
tmin<-ncvar_get(ftmin,"Band1")

bio<-numeric(2160*1080*19)
dim(bio)<-c(2160,1080,19)
for (i in 1:2160) {
  bio[i,,]<-biovars(prec[i,,],tmin[i,,],tmax[i,,])
}

Longvector<-ncvar_get(fprec,"lon")
Latvector<-ncvar_get(fprec,"lat")

dimX <- ncdim_def("lon", "degrees", Longvector)
dimY <- ncdim_def("lat", "degrees", Latvector)
dimT <- ncdim_def("bio", "unit", 1:19)
#    dimT <- dim.def.ncdf("Time", "days", 0, unlim = FALSE)

mv <- -9999 
bioclim <- ncvar_def( "bioclim", "unit",list(dimX,dimY,dimT),mv,prec="double")

nc = nc_create( "bioclim_plio.nc", bioclim )
#    nc <- create.ncdf( "writevals.nc", list(var1d,var2d,var3d))
ncvar_put(nc, bioclim, bio)

nc_close(nc)


    
