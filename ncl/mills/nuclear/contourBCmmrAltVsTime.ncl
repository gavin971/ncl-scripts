; This script creates Figure 1 for 2013 paper "Multi-decadal global cooling 
; and unprecedented ozone loss following a regional nuclear conflict" by 
; Mills, Toon, Lee-Taylor, and Robock. 

; Figure 1 caption: The time evolution of BC mass mixing ratio (kg BC/109 
; kg air) is shown for the average of the 5 Tg BC ensemble. The horizontal 
; axis shows time in years since the BC emission on January 1.

; Script written by Mike Mills
; Modified September 24, 2013

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

ostype = systemfunc("echo $OSTYPE")
if (ostype.eq."darwin") then
  path="$MODEL_PATH/cesm/runs/cesm1_0/"
else
  path="$MODEL_PATH/"
end if

exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"
exptFilepath=path+exptCase+"/h0gm/"+exptCase+".cam2.h0.globavg.nc"
exptFile=addfile(exptFilepath,"r")

BCmmr = exptFile->SOOT01
lev  = exptFile->lev
lat  = exptFile->lat
lon  = exptFile->lon
BCmmr!1="press"
BCmmr&press@units="hPa"
BCmmr&press@long_name="Pressure"
BCmmr&time=BCmmr&time/365.
BCmmr&time@units="Years"
BCmmr=BCmmr*1e9
BCmmr@units="ppbm"

;************************************************
; create plot
;************************************************
wks = gsn_open_wks("eps", "$IMG_PATH/nuclear/Fig01nuclearBCmmrAltVsTime")                ; open an eps file
;wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearBCmmrAltVsTime")                ; open a pdf file

; Change color map using named colors
colors = (/"white","black","white","cadetblue1","deepskyblue1","deepskyblue3","green4","greenyellow", \
           "yellow", "orange","orangered","red"/)
gsn_define_colormap(wks,colors)                  ; change colormap 
;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;gsn_define_colormap(wks,"BlAqGrYeOrRe")
 
res                       = True     ; plot mods desired
res@vpHeightF             = 0.6        ; change aspect ratio
;res@vpWidthF              = 1.0        ; change aspect ratio
  res@vpWidthF            = 1.0           ; change aspect ratio of plot
;  res@vpHeightF           = 0.3
;res@gsnDraw = False
;res@gsnFrame = False
res@gsnPaperOrientation = "landscape"
res@gsnMaximize         = True          ; expand plot size to maximum
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = True     ; turn on contour lines
res@cnLineLabelsOn        = True     ; turn on contour line labels
res@cnLineLabelFontHeightF=0.01
res@cnLevelFlag="LineAndLabel"
res@tmYLMode    = "Explicit" 
res@tmYLValues  = (/    1E-5 ,   1E-4  , 0.001 , 0.01 , 0.1 , 1., 10., 100./)
res@tmYLLabels  = (/"10~S~-5","10~S~-4","0.001","0.01","0.1","1","10","100"/)
res@tmXBMode    = "Manual" 
res@tmXBTickSpacingF  = 1.0
res@tmXBTickStartF = 1.0
res@cnLevelSelectionMode = "ExplicitLevels" 
res@cnLevels  = (/1,3,5,7,10,20,30,40,50/) 
res@trYMinF = 0.001
res@trXMaxF = 14.0
res@tmXBPrecision = 2
res@tmXBAutoPrecision = False
res@tmXBMinorOn = False
;res@trYMaxF = 500.
res@gsnLeftString=" "
res@gsnRightString=" "
res@tiXAxisString       = ""
res@tiXAxisString = "Time (years since BC injection)"
res@tiYAxisString = "Pressure (hPa)"
res@tiXAxisFontHeightF = 0.022
res@tiYAxisFontHeightF = 0.022
res@tiXAxisOffsetYF = 0.14
res@tiYAxisOffsetXF = 0.02
res@lbLabelBarOn         = True
;res@lbOrientation        = "Vertical"
res@lbOrientation        = "Horizontal"
res@lbTopMarginF         =  -3.8
res@lbBottomMarginF      =  3.9
res@lbRightMarginF       = -0.2
res@lbLeftMarginF        = 2.0
res@lbTitleOn = True
res@lbTitlePosition = "Bottom"
res@lbTitleString = "(kg BC/10~S~9~N~ kg air)"
res@lbTitleFontHeightF = 0.022
res@lbTitleOffsetF = -0.2
plot = gsn_csm_pres_hgt(wks,BCmmr(press|:, time|:),res)      ; contour the variable

delete(wks)  
