load "/usr/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/usr/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/usr/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/usr/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/usr/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "/usr/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  indir_prefix="/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/Indre/Esther_redescription_paleo/NorESM_Zhang/"
  indir_map="/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/mapping_data/"

  wks_type                         ="pdf"
  wks_type@wkPaperSize             ="A4"
;  wks_type@wkOrientation          ="landscape"
  wks                              = gsn_open_wks(wks_type, "test2")
    gsn_define_colormap(wks,"hotcolr_19lev")         ; choose color map
    i = NhlNewColor(wks,0.8,0.8,0.8)                   ; add gray to colormap

create_paleodata=False
write_site      =True

; 10ma2co2: 560 ppm CO2
; 10ma0608  350 ppm
; 20ma2co2  560 ppm CO2
; 20ma0329  420 ppm 
; 50ma2co2  560 ppm CO2
; 50ma1028  1120 ppm

jobname=(/"10ma2co2","10ma0608","20ma2co2","20ma0329","50ma2co2","50ma1028"/)
jobname2=(/"LateMio_Nor_560","LateMio_Nor_350","EarlyMio_Nor_560","EarlyMio_Nor_420","50ma2co2","50ma1028"/)

do m=0,3,1

if (create_paleodata) then

;### read paleodata

    fprec_pli           = addfile(indir_prefix+jobname(m)+"-PandT.nc","r")   
    ftas_pli            = addfile(indir_prefix+jobname(m)+"-PandT.nc","r")   
    ftmax_pli           = addfile(indir_prefix+jobname(m)+"-PandT.nc","r")   
    ftmin_pli           = addfile(indir_prefix+jobname(m)+"-PandT.nc","r") 
    fprec_mod           = addfile(indir_prefix+"PI-PandT.nc","r")   
    ftas_mod            = addfile(indir_prefix+"PI-PandT.nc","r")   
    ftmax_mod           = addfile(indir_prefix+"PI-PandT.nc","r")   
    ftmin_mod           = addfile(indir_prefix+"PI-PandT.nc","r") 

    pr_pli          = fprec_pli->PRECT(:,:,:)
    tas_pli         = fprec_pli->TREFHT(:,:,:)  
    tma_pli         = fprec_pli->TREFMXAV(:,:,:)
    tmi_pli         = fprec_pli->TREFMNAV(:,:,:)
    printVarSummary(tmi_pli)

    pr_mod          = fprec_mod->PRECT(:,:,:)
    tas_mod         = fprec_mod->TREFHT(:,:,:)    
    tma_mod         = fprec_mod->TREFMXAV(:,:,:)
    tmi_mod         = fprec_mod->TREFMNAV(:,:,:)
    printVarSummary(pr_mod)

;###calculate the delta-change of pliocene

    dpre_pli_rel  = (pr_pli-pr_mod)/ where(pr_mod.ne.0, pr_mod, default_fillvalue("float")) 
    dpre_pli_abs  = pr_pli-pr_mod
    dtas_pli = tas_pli - tas_mod
    dtma_pli = tma_pli - tma_mod
    dtmi_pli = tmi_pli - tmi_mod
    printVarSummary(dpre_pli_rel)
    printVarSummary(dpre_pli_abs)

;###regrid delta-change onto worldclim 10min grids
;###mapping files are generated on SAGA using ESMF:
;###Notes (02.2022): abel_backup is now kept in nird:/projects/NS2806K/huit 
;### To retrieve:  tar -xvf abel_backup.tar abel_backup/prepare_clm_input/ScriptGrid/SCRIPgrid_worldclim2_mask.nc
;                  scp abel_backup/prepare_clm_input/ScriptGrid/SCRIPgrid_worldclim2_mask.nc huit@saga.sigma2.no:/cluster/home/huit/
;                  scp SCRIPgrid_NorESM_nomask.nc huit@saga.sigma2.no:/cluster/home/huit/
;                  script grid of mode is kept in: /home/huitang/Documents/data/SCRIPgrid_data
;###/cluster/software/ESMF/8.0.0-intel-2019b/bin/ESMF_RegridWeightGen --ignore_unmapped -s SCRIPgrid_ECHAM5_nomask.nc -d /cluster/projects/nn2806k/huit/abel_backup/prepare_clm_input/ScriptGrid/SCRIPgrid_worldclim2_mask.nc -m bilinear -w map_echam5_to_1080x2160_nomask_blin_da_c200319.nc --64bit_offset

    Opt                = True
    Opt@PrintTimings   = True
    Opt@Debug          = True
    dpre_regrid_rel = ESMF_regrid_with_weights(dpre_pli_rel,indir_map+"/map_noresm_to_1080x2160_nomask_blin_da_c200319.nc",Opt)
    dpre_regrid_abs = ESMF_regrid_with_weights(dpre_pli_abs,indir_map+"/map_noresm_to_1080x2160_nomask_blin_da_c200319.nc",Opt)
    dtas_regrid = ESMF_regrid_with_weights(dtas_pli,indir_map+"/map_noresm_to_1080x2160_nomask_blin_da_c200319.nc",Opt)
    dtma_regrid = ESMF_regrid_with_weights(dtma_pli,indir_map+"/map_noresm_to_1080x2160_nomask_blin_da_c200319.nc",Opt)
    dtmi_regrid = ESMF_regrid_with_weights(dtmi_pli,indir_map+"/map_noresm_to_1080x2160_nomask_blin_da_c200319.nc",Opt)

    printVarSummary(dpre_regrid_rel)
    printVarSummary(dpre_regrid_abs)
    printVarSummary(dtas_regrid)

;### add delta-change to worldclim 10min data

    fprec           = addfile("/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_prec/wc2.0_10m_prec.nc","r")   
    ftas            = addfile("/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tavg/wc2.0_10m_tavg.nc","r")   
    ftmax           = addfile("/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tmax/wc2.0_10m_tmax.nc","r")   
    ftmin           = addfile("/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tmin/wc2.0_10m_tmin.nc","r")   

    lat             = fprec->lat(:)
    lon             = fprec->lon(:)
    printVarSummary(lat)
    printVarSummary(lon)

    pre_wc          = tofloat(fprec->Band1(:,:,:))
    tas_wc          = ftas->Band1(:,:,:)    
    tmax_wc         = ftmax->Band1(:,:,:)
    tmin_wc         = ftmin->Band1(:,:,:)
    printVarSummary(pre_wc)

;#there are many ways to add anomoly fro precipition, (1)absoluted difference (2) relative difference (3) mix of absolute and relative
;#here, I use combined approach, which has been used in my phd study; and also relative difference
    pre_wc_pli      = pre_wc + dpre_regrid_abs
    do t=0,11,1
       do i=0,2159,1
          do j=0,1079,1 
             if ((.not.ismissing(dpre_regrid_rel(t,j,i))) .and. (dpre_regrid_rel(t,j,i) .lt. 0)) then
                 pre_wc_pli(t,j,i)      = pre_wc(t,j,i) + dpre_regrid_rel(t,j,i) * pre_wc(t,j,i)
             end if
             if ( (.not.ismissing(pre_wc_pli(t,j,i))) .and. (pre_wc_pli(t,j,i)  .lt. 0)) then
                 pre_wc_pli(t,j,i)      = 0
             end if
          end do
       end do
    end do

    tas_wc_pli      = tas_wc + dtas_regrid
    tmax_wc_pli     = tmax_wc + dtma_regrid
    tmin_wc_pli     = tmin_wc + dtmi_regrid

    printVarSummary(pre_wc_pli)

;###Put the new data into netcdf file

   system("mkdir -p "+indir_prefix+jobname(m))
   system("cp /media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_prec/wc2.0_10m_prec.nc "+indir_prefix+jobname(m)+"/wc2.0_10m_prec_plio.nc")
   system("cp /media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tavg/wc2.0_10m_tavg.nc "+indir_prefix+jobname(m)+"/wc2.0_10m_tavg_plio.nc")
   system("cp /media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tmax/wc2.0_10m_tmax.nc "+indir_prefix+jobname(m)+"/wc2.0_10m_tmax_plio.nc")
   system("cp /media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/WorldClim/WorldClim2.0/wc2.0_10m_tmin/wc2.0_10m_tmin.nc "+indir_prefix+jobname(m)+"/wc2.0_10m_tmin_plio.nc")
   
   fp   = addfile(indir_prefix+jobname(m)+"/wc2.0_10m_prec_plio.nc","w")      
   fta  = addfile(indir_prefix+jobname(m)+"/wc2.0_10m_tavg_plio.nc","w")  
   ftma = addfile(indir_prefix+jobname(m)+"/wc2.0_10m_tmax_plio.nc","w")  
   ftmi = addfile(indir_prefix+jobname(m)+"/wc2.0_10m_tmin_plio.nc","w")  

   fp->Band1   = (/toshort(pre_wc_pli(:,:,:))/)  
   fta->Band1  = (/tas_wc_pli(:,:,:)/)  
   ftma->Band1 = (/tmax_wc_pli(:,:,:)/)  
   ftmi->Band1 = (/tmin_wc_pli(:,:,:)/)  
   
end if   


;###calculate bioclim variables with new data
;# See the R script: bioclim_var_echam.r


if (write_site) then
;###interpolate new bioclim variables onto fossil sites and write the output
    proxy        = asciiread("/media/huitang/f1e4880e-3c1f-4fea-8824-a6374ea26bf6/huitang/Documents/data/Indre/indre_scripts/Esther_paleo/fossil_locality/sites_china_updates_7",(/171,3/),"double")
    fbio         = addfile(indir_prefix+jobname(m)+"/bioclim_plio.nc","r") 
;#read localities
    proxynum       = proxy(:,0)
    proxylat       = proxy(:,1)     
    proxylon       = proxy(:,2)
    bioclim        = fbio->bioclim(:,:,:)
    lat            = fbio->lat(:)
    lon            = fbio->lon(:)

    ;printVarSummary(proxylat)
    print(proxylat)
    print(proxylon)

    numd           = dimsizes(proxylat)        ; number of distinct markers
    print(numd)

    var_proxy      = new((/numd,22/),"double") 

;#data interpolation

   jlat = ind_nearest_coord (proxylat,lat,0)
   ilon = ind_nearest_coord (proxylon,lon,0)
   printVarSummary(jlat)

   do q=0,numd-1,1
    var_proxy(q,0)=proxynum(q)
    var_proxy(q,1)=proxylon(q)
    var_proxy(q,2)=proxylat(q)
    do mon=0,18,1    
      var_proxy(q,3+mon)=(/bioclim(mon,jlat(q),ilon(q))/)
    end do
   end do
 
    var_proxy@_FillValue=-9999

    printVarSummary(var_proxy)

    fmtf="22F20.2"
    opt =True
    opt@fout="bioclim_"+jobname2(m)+".txt"
    opt@row =False
    
    write_matrix(var_proxy,fmtf,opt)

end if

end do

end
