load "/usr/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/usr/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/usr/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/usr/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/usr/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "/usr/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  indir_prefix="/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_bio"

  create_paleodata=False
  write_site      =True

  if (write_site) then
  ;###interpolate new bioclim variables onto fossil sites and write the output
    proxy        = asciiread("/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/Indre/indre_scripts/Esther_paleo/fossil_locality/sites_china_updates_7",(/171,3/),"double")
    fbio         = addfile(indir_prefix+"/bio01-19.nc","r") 
  ;#read localities
    proxynum       = proxy(:,0)
    proxylat       = proxy(:,1)     
    proxylon       = proxy(:,2)
    lat            = fbio->lat(:)
    lon            = fbio->lon(:)
    nlat           = dimsizes(lat)
    nlon           = dimsizes(lon)
    bioclim        = new((/19,nlat,nlon/),double)
    bioclim(0,:,:) = fbio->bio01(:,:)
    bioclim(1,:,:) = fbio->bio02(:,:)
    bioclim(2,:,:) = fbio->bio03(:,:)
    bioclim(3,:,:) = fbio->bio04(:,:)
    bioclim(4,:,:) = fbio->bio05(:,:)
    bioclim(5,:,:) = fbio->bio06(:,:)
    bioclim(6,:,:) = fbio->bio07(:,:)
    bioclim(7,:,:) = fbio->bio08(:,:)
    bioclim(8,:,:) = fbio->bio09(:,:)
    bioclim(9,:,:) = fbio->bio10(:,:)
    bioclim(10,:,:)= fbio->bio11(:,:)
    bioclim(11,:,:)= fbio->bio12(:,:)
    bioclim(12,:,:)= fbio->bio13(:,:)
    bioclim(13,:,:)= fbio->bio14(:,:)
    bioclim(14,:,:)= fbio->bio15(:,:)
    bioclim(15,:,:)= fbio->bio16(:,:)
    bioclim(16,:,:)= fbio->bio17(:,:)
    bioclim(17,:,:)= fbio->bio18(:,:)
    bioclim(18,:,:)= fbio->bio19(:,:)

    print(proxylat)
    print(proxylon)


    numd           = dimsizes(proxylat)        ; number of distinct markers
    print(numd)

    var_proxy      = new((/numd,22/),"double") 

   ;#data interpolation

   jlat = ind_nearest_coord (proxylat,lat,0)
   ilon = ind_nearest_coord (proxylon,lon,0)
   printVarSummary(jlat)

   do m=0,numd-1,1
    var_proxy(m,0)=proxynum(m)
    var_proxy(m,1)=proxylon(m)
    var_proxy(m,2)=proxylat(m)
    do mon=0,18,1    
      var_proxy(m,3+mon)=(/bioclim(mon,jlat(m),ilon(m))/)
    end do
   end do
 
    var_proxy@_FillValue=-9999

    printVarSummary(var_proxy)

    fmtf="22F20.2"
    opt =True
    opt@fout="bioclim_present.txt"
    opt@row =False
    
    write_matrix(var_proxy,fmtf,opt)

end if

end 
