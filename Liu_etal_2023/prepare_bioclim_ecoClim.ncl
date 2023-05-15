load "/usr/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/usr/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/usr/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/usr/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/usr/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "/usr/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  indir_prefix="/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/Paleo_bioclim/EcoClim"

  wks_type                         ="pdf"
  wks_type@wkPaperSize             ="A4"
;  wks_type@wkOrientation          ="landscape"
  wks                              = gsn_open_wks(wks_type, "test2")
    gsn_define_colormap(wks,"hotcolr_19lev")         ; choose color map
    i = NhlNewColor(wks,0.8,0.8,0.8)                   ; add gray to colormap

create_paleodata=False
write_site      =True


if (create_paleodata) then
;### read paleodata

    pr_pli   = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/pr_CCSM_Plio3Ma.txt",(/250560,14/),"float")
    tas_pli        = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/tas_CCSM_Plio3Ma.txt",(/250560,14/),"float")
    tma_pli        = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/tasmax_CCSM_Plio3Ma.txt",(/250560,14/),"float")
    tmi_pli        = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/tasmin_CCSM_Plio3Ma.txt",(/250560,14/),"float")
    pr_mod   = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/pr_CCSM_Modern.txt",(/250560,14/),"float")
    tas_mod       = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/tas_CCSM_Modern.txt",(/250560,14/),"float")
    tma_mod        = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/tasmax_CCSM_Modern.txt",(/250560,14/),"float")
    tmi_mod        = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/tasmin_CCSM_Modern.txt",(/250560,14/),"float")
    
    printVarSummary(pr_pli(:,12))   
;    print(tas_pli(:,0))   
;    print(tma_pli(:,0))   
;    print(tmi_pli(:,0))   
;    print(pr_mod(:,13))   
;    print(tas_mod(:,13))   
;   print(tma_mod(:,13))   
;    print(tmi_mod(:,13))    

;###calculate the delta-change of pliocene

    d_pr_pli_rel  = (pr_pli(:,2:13)-pr_mod(:,2:13))/pr_mod(:,2:13)
    d_pr_pli_abs  = pr_pli(:,2:13)-pr_mod(:,2:13)
    d_tas_pli = tas_pli(:,2:13) - tas_mod(:,2:13)
    d_tma_pli = tma_pli(:,2:13) - tma_mod(:,2:13)
    d_tmi_pli = tmi_pli(:,2:13) - tmi_mod(:,2:13)

    lon_pli   = onedtond(pr_pli(:,0),(/720,348/))
    lat_pli   = onedtond(pr_pli(:,1),(/720,348/))
    print(lat_pli(0,:))
    print(lon_pli(:,0))
    lon_pli1  = lon_pli(:,0)
    lat_pli1  = lat_pli(0,:)       

    dpre_pli_rel=new((/12,348,720/),"float")
    dpre_pli_abs=new((/12,348,720/),"float")
    dtas_pli=new((/12,348,720/),"float")
    dtma_pli=new((/12,348,720/),"float")
    dtmi_pli=new((/12,348,720/),"float")
 
    do t=0,11,1
       do i=0,719,1
          do j=0,347,1          
             dpre_pli_rel(t,j,i)= d_pr_pli_rel(i*348+j,t)
             dpre_pli_abs(t,j,i)= d_pr_pli_abs(i*348+j,t)
             dtas_pli(t,j,i)= d_tas_pli(i*348+j,t)
             dtma_pli(t,j,i)= d_tma_pli(i*348+j,t)
             dtmi_pli(t,j,i)= d_tmi_pli(i*348+j,t)
          end do
       end do
    end do


;###regrid delta-change onto worldclim 10min grids
;###mapping files are generated on SAGA using ESMF:
;###/cluster/software/ESMF/8.0.0-intel-2019b/bin/ESMF_RegridWeightGen --ignore_unmapped -s SCRIPgrid_ecoClim_nomask_rect.nc -d /cluster/projects/nn2806k/huit/abel_backup/prepare_clm_input/ScriptGrid/SCRIPgrid_worldclim2_mask.nc -m bilinear -w map_ecoclimrec_to_1080x2160_nomask_blin_da_c200319.nc --64bit_offset
;###/cluster/software/ESMF/8.0.0-intel-2019b/bin/ESMF_RegridWeightGen --ignore_unmapped -s SCRIPgrid_ecoClim_nomask.nc -d /cluster/projects/nn2806k/huit/abel_backup/prepare_clm_input/ScriptGrid/SCRIPgrid_worldclim2_mask.nc -m bilinear -w map_ecoclim_to_1080x2160_nomask_blin_da_c200319.nc --64bit_offset --src_loc center (This does not work, too slow)

    Opt                = True
    Opt@PrintTimings   = True
    Opt@Debug          = True
    dpre_regrid_rel = ESMF_regrid_with_weights(dpre_pli_rel,"/home/huitang/Documents/data/Paleo_bioclim/map_ecoclimrec_to_1080x2160_nomask_blin_da_c200319.nc",Opt)
    dpre_regrid_abs = ESMF_regrid_with_weights(dpre_pli_abs,"/home/huitang/Documents/data/Paleo_bioclim/map_ecoclimrec_to_1080x2160_nomask_blin_da_c200319.nc",Opt)
    dtas_regrid = ESMF_regrid_with_weights(dtas_pli,"/home/huitang/Documents/data/Paleo_bioclim/map_ecoclimrec_to_1080x2160_nomask_blin_da_c200319.nc",Opt)
    dtma_regrid = ESMF_regrid_with_weights(dtma_pli,"/home/huitang/Documents/data/Paleo_bioclim/map_ecoclimrec_to_1080x2160_nomask_blin_da_c200319.nc",Opt)
    dtmi_regrid = ESMF_regrid_with_weights(dtmi_pli,"/home/huitang/Documents/data/Paleo_bioclim/map_ecoclimrec_to_1080x2160_nomask_blin_da_c200319.nc",Opt)

    printVarSummary(dpre_regrid_rel)
    printVarSummary(dtas_regrid)

;### add delta-change to worldclim 10min data

    fprec           = addfile("/home/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_prec/wc2.0_10m_prec.nc","r")   
    ftas            = addfile("/home/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tavg/wc2.0_10m_tavg.nc","r")   
    ftmax           = addfile("/home/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tmax/wc2.0_10m_tmax.nc","r")   
    ftmin           = addfile("/home/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tmin/wc2.0_10m_tmin.nc","r")   

    lat             = fprec->lat(:)
    lon             = fprec->lon(:)
    printVarSummary(lat)
    printVarSummary(lon)

    pre_wc          = tofloat(fprec->Band1(:,:,:))
    tas_wc          = ftas->Band1(:,:,:)    
    tmax_wc         = ftmax->Band1(:,:,:)
    tmin_wc         = ftmin->Band1(:,:,:)

;#there are many ways to add anomoly fro precipition, (1)absoluted difference (2) relative difference (3) mix of absolute and relative
;#here, I use combined approach, which has been used in my phd study; and also relative difference
    pre_wc_pli      = pre_wc + dpre_regrid_abs
    do t=0,11,1
       do i=0,2159,1
          do j=0,1079,1 
             if ((.not.ismissing(dpre_regrid_abs(t,j,i))) .and. (dpre_regrid_abs(t,j,i) .ge. 0)) then
                 pre_wc_pli(t,j,i)      = pre_wc(t,j,i) + dpre_regrid_abs(t,j,i)
             else
                 pre_wc_pli(t,j,i)      = pre_wc(t,j,i) + dpre_regrid_rel(t,j,i) * pre_wc(t,j,i)
             end if
          end do
       end do
    end do

    tas_wc_pli      = tas_wc + dtas_regrid
    tmax_wc_pli     = tmax_wc + dtma_regrid
    tmin_wc_pli     = tmin_wc + dtmi_regrid

;###Put the new data into netcdf file

   system("cp /home/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_prec/wc2.0_10m_prec.nc /home/huitang/Documents/data/Paleo_bioclim/EcoClim/wc2.0_10m_prec_ecoclim_plio.nc")
   system("cp /home/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tavg/wc2.0_10m_tavg.nc /home/huitang/Documents/data/Paleo_bioclim/EcoClim/wc2.0_10m_tavg_ecoclim_plio.nc")
   system("cp /home/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tmax/wc2.0_10m_tmax.nc /home/huitang/Documents/data/Paleo_bioclim/EcoClim/wc2.0_10m_tmax_ecoclim_plio.nc")
   system("cp /home/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tmin/wc2.0_10m_tmin.nc /home/huitang/Documents/data/Paleo_bioclim/EcoClim/wc2.0_10m_tmin_ecoclim_plio.nc")
   fp   = addfile("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/wc2.0_10m_prec_ecoclim_plio.nc","w")      
   fta  = addfile("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/wc2.0_10m_tavg_ecoclim_plio.nc","w")  
   ftma = addfile("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/wc2.0_10m_tmax_ecoclim_plio.nc","w")  
   ftmi = addfile("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/wc2.0_10m_tmin_ecoclim_plio.nc","w")  

   fp->Band1   = (/toshort(pre_wc_pli(:,:,:))/)  
   fta->Band1  = (/tas_wc_pli(:,:,:)/)  
   ftma->Band1 = (/tmax_wc_pli(:,:,:)/)  
   ftmi->Band1 = (/tmin_wc_pli(:,:,:)/)  
   
end if   
;###calculate bioclim variables with new data

;# See the R script: bioclim_var_ecoclim.r


if (write_site) then
;###interpolate new bioclim variables onto fossil sites and write the output
    proxy        = asciiread("/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/Indre/indre_scripts/Esther_paleo/fossil_locality/sites_china_updates_7",(/171,3/),"double")
    fbio          = addfile("/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/Paleo_bioclim/EcoClim/bioclim_ecoclim_plio.nc","r") 
;#read localities
    proxynum       = proxy(:,0)
    proxylat       = proxy(:,1)     
    proxylon       = proxy(:,2)
    bioclim        = fbio->bioclim(:,:,:)
    lat            = fbio->lat(:)
    lon            = fbio->lon(:)

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
    opt@fout="bioclim_MidPlio_405.txt"
    opt@row =False
    
    write_matrix(var_proxy,fmtf,opt)

end if
    
end 
