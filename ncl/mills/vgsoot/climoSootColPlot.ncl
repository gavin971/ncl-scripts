load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotClimoSootCol")
procedure plotClimoSootCol(season:string, rootPath:string) 
local wks,res,lev
begin 

  fileName=rootPath+"vgsoot08_"+season+"_climo.nc"
  a = addfile(fileName,"r")
  SOOTMD_col=a->SOOTMD_col

  lat=a->lat
  lon=a->lon
  
  xdim=dimsizes(SOOTMD_col)
  print("xdim="+xdim)
  ntimes=min((/xdim(0)/))
    
  wks = gsn_open_wks("pdf", "$IMAGE_PATH/ClimoSootColumnPolarMap"+season)
;  gsn_define_colormap(wks,"BlRe")
;  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  gsn_define_colormap(wks,"rainbow")
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
;  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "portrait"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7

  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=500
  res@cnMaxLevelValF=9000
  res@cnLevelSpacingF     = 500.
  
  res@gsnCenterString      = "Soot Column (ng m~S~-2~N~)"
;  res@tiMainString         =   
  res@gsnRightString      = season
  res@lbLabelAngleF       = 60
  res@gsnPolar   = "NH"                          ; specify the hemisphere

; Create and draw a contour plot over a polar stereographic map.
  plot = gsn_csm_contour_map_polar(wks,SOOTMD_col(0,:,:),res)   ; plot model data 
   
  res@gsnPolar   = "SH"                          ; specify the hemisphere

; Create and draw a contour plot over a polar stereographic map.
  plot = gsn_csm_contour_map_polar(wks,SOOTMD_col(0,:,:),res)   ; plot model data  
  
end


;season="ANN"
;season="DJF"
season="JJA"

rootPath="/Volumes/Data/Models/waccm/run/carma-soot/vgsoot08/climo/2035-2064/"

plotClimoSootCol(season,rootPath)
