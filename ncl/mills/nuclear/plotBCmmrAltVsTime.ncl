load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

path="$MODEL_PATH/ccsm/case/"

exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"
cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"

exptFilepath=path+exptCase+"/h0gm/"+exptCase+".cam2.h0.globavg.nc"
cntlFilepath=path+cntlCase+"/h0gm/"+cntlCase+".cam2.h0.globavg.nc"

exptFile=addfile(exptFilepath,"r")
cntlFile=addfile(cntlFilepath,"r")

BCmmr = exptFile->SOOT01
Te    = exptFile->T
Tc    = cntlFile->T
Tdiff = Te ; get metadata
Tdiff = Tdiff - Tc
;ps   = exptFile->PS
;p0   = exptFile->P0
;hyam = exptFile->hyam
;hybm = exptFile->hybm
;time = exptFile->time
lev  = exptFile->lev
lat  = exptFile->lat
lon  = exptFile->lon
;pm   = (ps,p0,hyai,hybi)
BCmmr!1="press"
BCmmr&press@units="hPa"
BCmmr&press@long_name="Pressure"
BCmmr&time=BCmmr&time/365.
BCmmr&time@units="Years"
BCmmr=BCmmr*1e9
BCmmr@units="ppbm"

Tdiff!1="press"
Tdiff&press=BCmmr&press
Tdiff&time=BCmmr&time

;************************************************
; create plot
;************************************************
;wks = gsn_open_wks("eps", "$IMG_PATH/nuclearBCmmrAltVsTime")                ; open an eps file
wks = gsn_open_wks("pdf", "$IMG_PATH/nuclearBCmmrAltVsTime")                ; open a pdf file
;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
gsn_define_colormap(wks,"hotres")

nPanel=2
plot = new(nPanel, graphic)

res                       = True     ; plot mods desired
res@vpHeightF             = 0.5 / nPanel              ; change aspect ratio
res@gsnDraw = False
res@gsnFrame = False
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
res@cnLevels  = (/1,2,3,5,7,10,20,30,50,70,100/) 
res@trYMinF = 0.001
res@trYMaxF = 500.
res@gsnLeftString=" "
res@gsnRightString=" "
res@tiXAxisString       = ""
res@lbLabelBarOn         = True
res@lbOrientation        = "Vertical"
res@lbTopMarginF         =  0.0
res@lbBottomMarginF      = -3.0
res@lbRightMarginF       = -0.1
res@lbLeftMarginF        = 0.3
plot(0) = gsn_csm_pres_hgt(wks,BCmmr(press|:, time|:),res)      ; contour the variable

res@lbLabelBarOn         = True
res@lbTopMarginF         = 5
res@lbBottomMarginF      = -7
res@tiXAxisString       = "Time (years)"
plot(1) = gsn_csm_pres_hgt(wks,Tdiff(press|:, time|:),res)      ; contour the variable

resP = True
resP@gsnMaximize         = True          ; expand plot size to maximum
resP@gsnPaperOrientation = "landscape"
gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
delete(plot)  
