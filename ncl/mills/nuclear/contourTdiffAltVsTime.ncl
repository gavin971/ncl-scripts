; This script creates Figure 5 for 2013 paper "Multi-decadal global cooling 
; and unprecedented ozone loss following a regional nuclear conflict" by 
; Mills, Toon, Lee-Taylor, and Robock. 

; Figure 5 caption:  

; Script written by Mike Mills
; Modified September 24, 2013

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

path="$MODEL_PATH/cesm/runs/cesm1_0/"

exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"
cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"

exptFilepath=path+exptCase+"/h0gm/"+exptCase+".cam2.h0.globavg.nc"
cntlFilepath=path+cntlCase+"/h0gm/"+cntlCase+".cam2.h0.globavg.nc"

exptFile=addfile(exptFilepath,"r")
cntlFile=addfile(cntlFilepath,"r")

Te    = exptFile->T
Tc    = cntlFile->T
Tdiff = Te ; get metadata
Tdiff = Tdiff - Tc
lev  = exptFile->lev
lat  = exptFile->lat
lon  = exptFile->lon

Tdiff!1="press"
Tdiff&press@units="hPa"
Tdiff&press@long_name="Pressure"
Tdiff&time=Tdiff&time/365.
Tdiff&time@units="Years"

;************************************************
; create plot
;************************************************
wks = gsn_open_wks("eps", "$IMG_PATH/nuclear/Fig05_contourTdiffAltVsTime") ; open a pdf file
;wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearTdiffAltVsTime")             ; open a pdf file

colors = (/"white","black","lightblue","white","yellow1", \
           "yellow3","orange1","orange2","orangered","red1","red2","red4"/)
gsn_define_colormap(wks,colors)                  ; change colormap 

res                       = True     ; plot mods desired
res@vpHeightF             = 0.6              ; change aspect ratio
res@vpWidthF            = 1.0           ; change aspect ratio of plot
res@gsnPaperOrientation = "landscape"
res@gsnMaximize         = True          ; expand plot size to maximum
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = True     ; turn on contour lines
res@cnLineLabelsOn        = True     ; turn on contour line labels
res@cnLineLabelFontHeightF=0.01
res@tmYLMode    = "Explicit" 
res@tmYLValues  = (/    1E-5 ,   1E-4  , 0.001 , 0.01 , 0.1 , 1., 10., 100./)
res@tmYLLabels  = (/"10~S~-5","10~S~-4","0.001","0.01","0.1","1","10","100"/)
res@tmXBMode    = "Manual" 
res@tmXBTickSpacingF  = 1.0
res@tmXBTickStartF = 1.0
res@cnLevelSelectionMode = "ExplicitLevels" 
res@cnLevels  = (/-1,1,3,5,7,10,30,50,70/) 
res@trYMinF = 0.001
res@trXMaxF = 20.0
res@tmXBPrecision = 2
res@tmXBAutoPrecision = False
res@tmXBMinorOn = False
;res@trYMaxF = 500.
res@gsnLeftString=" "
res@gsnRightString=" "
res@tiXAxisString       = ""
res@lbLabelBarOn         = True
res@lbOrientation        = "Vertical"
res@lbTopMarginF         = -0.2
res@lbBottomMarginF      = -0.5
res@lbRightMarginF       = 0.3
res@lbLeftMarginF        = 0.1

;res@lbTopMarginF         = 5
;res@lbBottomMarginF      = -7
res@tiXAxisString       = "Time (years since BC injection)"
res@tiYAxisString       = "Pressure (hPa)"
res@tiXAxisFontHeightF = 0.022
res@tiYAxisFontHeightF = 0.022
res@tiYAxisOffsetXF = 0.02
plot = gsn_csm_pres_hgt(wks,Tdiff(press|:, time|:),res)      ; contour the variable

delete(wks)
