load "/usr/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/usr/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "/usr/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/usr/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/usr/lib/ncarg/nclscripts/csm/popRemap.ncl"
load "/usr/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  indir_prefix="/home/huitang/Documents/data/Indre/Esther_redescription_paleo/00360ctrl_echam5"

  wks_type                         ="pdf"
  wks_type@wkPaperSize             ="A4"
;  wks_type@wkOrientation          ="landscape"
  wks                              = gsn_open_wks(wks_type, "test2")
    gsn_define_colormap(wks,"hotcolr_19lev")         ; choose color map
    i = NhlNewColor(wks,0.8,0.8,0.8)                   ; add gray to colormap


;**********
;draw figures
;**********
   f1           = addfile(indir_prefix+"/Mean_Monthly_Sums_precip_tot_2492-2500.nc","r")
    a           = addfile("/usr/lib/ncarg/data/cdf/landsea.nc","r")
;   proxy        = asciiread(indir_prefix+"/global100kmsquares.csv",(/51468,3/),"float")  ; 

;----read era data
   lat          = f1->lat(:)
   lon          = f1->lon(:)

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

   rectilinear_to_SCRIP(indir_prefix+"/SCRIPgrid_ECHAM5_nomask.nc",lat,lon,Opt1)
    
end 
