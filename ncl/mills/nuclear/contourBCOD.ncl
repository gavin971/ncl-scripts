; This script creates Figure 2 for 2013 paper "Multi-decadal global cooling 
; and unprecedented ozone loss following a regional nuclear conflict" by 
; Mills, Toon, Lee-Taylor, and Robock. 

; Figure 2 caption: The time evolution of zonal mean BC optical depth is shown 
; for the 5 Tg BC ensemble average. The vertical axis shows latitude. The 
; horizontal average shows time in years.

; Script written by Mike Mills
; Modified September 24, 2013

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

ostype = systemfunc("echo $OSTYPE")
if (ostype.eq."darwin") then
  path="$MODEL_PATH/cesm/runs/cesm1_0/"
else
  path="$MODEL_PATH/"
end if

exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"
exptFilepath=path+exptCase+"/proc/"+exptCase+".h0zm.ODV_SOOT01.nc"
exptFile=addfile(exptFilepath,"r")

BCOD=exptFile->ODV_SOOT01
BCOD&time=BCOD&time/365.
BCOD&time@units="Years"
  
  
  pdffile = "$IMAGE_PATH/nuclear/Fig02ContourBCOD."+exptCase
  wks = gsn_open_wks("eps", pdffile)
;  wks = gsn_open_wks("pdf", pdffile)
  
;  gsn_define_colormap(wks,"posneg_1")
;gsn_define_colormap(wks,"BlWhRe")                  ; change colormap 
; Change color map using named colors
  colors = (/"white","black","white","yellow1","yellow2","yellow3", \
             "orange1","orange2","orange3","orangered","red1","red3","red4" /)
  gsn_define_colormap(wks,colors)                  ; change colormap 

  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
  res@cnLinesOn             = True     ; turn on contour lines
  res@cnLineLabelsOn        = True     ; turn on contour line labels
  res@cnLineLabelFontHeightF=0.012
;  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = False          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.55
;  res@gsnLeftStringFontHeightF = 24
;  res@gsnRightStringFontHeightF = 24
  res@tiXAxisString       = "Time (years since BC injection)"
;  res@tiXAxisOffsetYF = 0.12
  res@tiYAxisString       = "Latitude"
  res@tiXAxisFontHeightF = 0.025
  res@tiYAxisFontHeightF = 0.025
  res@cnLineLabelPlacementMode= "Computed"
  res@cnLineLabelDensityF = 5
  res@tmXBLabelFontHeightF =0.025
  res@tmYLLabelFontHeightF =0.025
 
;  res@tmXBMaxTicks      = 10
;  res@tmXBMinorPerMajor = 3
;  res@tmXBAutoPrecision = False
;  res@tmXBPrecision = 2
  res@tmXBMode    = "Manual" 
  res@tmXBTickSpacingF  = 1.0
  res@tmXBTickStartF = 1.0

res@trXMaxF = 16.0
res@tmXBPrecision = 2
res@tmXBAutoPrecision = False
res@tmXBMinorOn = False

;===================================================================  
    
; Define our own contour levels.  
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=120
  res@cnMaxLevelValF=510
  res@cnLevelSpacingF     = 30.
;  res@lbLabelAngleF       = 60
  res@gsnLeftString       = "Control run"
;  plot = gsn_csm_lat_time(wks,O3_col_c(lat|:, time|:),res)      ; contour the variable
  res@gsnLeftString       = "Soot run"

  res@cnLevelSelectionMode = "ExplicitLevels" 
;  delete(res@cnLevels)
  res@cnLevels          = (/0.01,0.02,0.03,0.04,0.05,0.07,0.09,0.11,0.13,0.15/)
;  delete(res@cnLineThicknesses)
;  res@cnLineThicknesses = (/ 1,  1 , 1 , 1 , 1 , 1 , 1 , 1, 1,4, 1, 1/)
  res@cnMonoLineThickness = False
  res@gsnLeftString       = ""
  res@gsnRightString      = ""
  res@gsnCenterString      = "BC Optical Depth"
  res@gsnCenterStringFontHeightF   = 0.025
  res@lbLabelFontHeightF   = 0.025
  res@lbOrientation        = "Vertical"
  res@lbTopMarginF         =  -0.2
  res@lbBottomMarginF      =  -0.6
  res@lbRightMarginF       = 0.6
  res@lbLeftMarginF        = -0.05
  res@tmYLLabels=(/"90~S~o~N~S","60~S~o~N~S","30~S~o~N~S","0","30~S~o~N~N",\
    "60~S~o~N~N","90~S~o~N~N"/)
  plot = gsn_csm_lat_time(wks,BCOD(lat|:, time|:),res)      ; contour the variable
  delete(plot)
  
  ostype = systemfunc("echo $OSTYPE")
  if (ostype.eq."darwin") then
    opencmd="open "+pdffile+".pdf"
    print(opencmd)
    open=systemfunc(opencmd)
    print(open)
  end if
;end
