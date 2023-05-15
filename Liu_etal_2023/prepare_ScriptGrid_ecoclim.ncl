load "/usr/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/usr/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/usr/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/usr/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/usr/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "/usr/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  indir_prefix="/home/huitang/Documents/data/Paleo_bioclim/EcoClim"

  wks_type                         ="pdf"
  wks_type@wkPaperSize             ="A4"
;  wks_type@wkOrientation          ="landscape"
  wks                              = gsn_open_wks(wks_type, "test2")
    gsn_define_colormap(wks,"hotcolr_19lev")         ; choose color map
    i = NhlNewColor(wks,0.8,0.8,0.8)                   ; add gray to colormap


;**********
;draw figures
;**********
    pr_pli   = asciiread("/home/huitang/Documents/data/Paleo_bioclim/EcoClim/pr_CCSM_Plio3Ma.txt",(/250560,14/),"float")
    a           = addfile("/usr/lib/ncarg/data/cdf/landsea.nc","r")

;----read era data
    lon_pli   = onedtond(pr_pli(:,0),(/720,348/))
    lat_pli   = onedtond(pr_pli(:,1),(/720,348/))
    print(lat_pli(0,:))
    print(lon_pli(:,0))
    lon       = lon_pli(:,0)
    lat       = lat_pli(0,:)    

    lsdata = a->LSMASK
    lsm  = tointeger(landsea_mask(lsdata,lat,lon))
    printVarSummary(lsm)
    lsm  = where (lsm.eq.2, 1, lsm) ; make lakes "land"
    lsm  = where (lsm.eq.3, 1, lsm) ; make island "land"
    lsm  = where (lsm.eq.4, 1, lsm) ; make ice "land"

   printVarSummary(lat)
   printVarSummary(lon)
   printVarSummary(lsm)
 
   Opt1                = True
   Opt1@ForceOverwrite = True
;   Opt1@GridMask       = lsm
;    Opt1@GridMask         = (/1/)
;    Opt1@GridCornerLat  = (/yfloor,yfloor,yceil,yceil/)
;    Opt1@GridCornerLon  = (/xfloor,xceil,xceil,xfloor/)

   rectilinear_to_SCRIP(indir_prefix+"/SCRIPgrid_ecoClim_nomask_rect.nc",lat,lon,Opt1)
 
;    unstructured_to_ESMF(indir_prefix+"/SCRIPgrid_ecoClim_nomask.nc",pr_pli(:,1),pr_pli(:,0),Opt1)   
end 
