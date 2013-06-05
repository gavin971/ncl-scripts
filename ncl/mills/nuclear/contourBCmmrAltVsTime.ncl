load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

path="$MODEL_PATH/"

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
;wks = gsn_open_wks("eps", "$IMG_PATH/nuclearBCmmrAltVsTime")                ; open an eps file
wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearBCmmrAltVsTime")                ; open a pdf file

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
;res@tiXAxisString = "Time (years)"
res@lbLabelBarOn         = True
res@lbOrientation        = "Horizontal"
;res@lbTopMarginF         =  0.0
;res@lbBottomMarginF      =  0.0
;res@lbRightMarginF       = -0.1
;res@lbLeftMarginF        = 0.3
plot = gsn_csm_pres_hgt(wks,BCmmr(press|:, time|:),res)      ; contour the variable

delete(wks)  
