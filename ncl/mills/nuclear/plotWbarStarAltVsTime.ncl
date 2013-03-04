load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

path="$MODEL_PATH/ccsm/case/"

exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"
cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"

exptFilepath=path+exptCase+"/proc/concat_tem_diag0_"+exptCase+".nc"
cntlFilepath=path+cntlCase+"/proc/concat_tem_diag0_"+cntlCase+".nc"

exptFile=addfile(exptFilepath,"r")
cntlFile=addfile(cntlFilepath,"r")

eWstar = exptFile->WSTAR
cWstar = cntlFile->WSTAR

time   = exptFile->time
time = time/365.
time@units="Years"
press = exptFile->lev
press@units="hPa"
press@long_name="Pressure"

eWstar!1="press"
eWstar&press=press
eWstar&time=time

cWstar!1="press"
cWstar&press=press
cWstar&time=time

eWstar_tropics=dim_avg_Wrap(eWstar(:,:,{-22:22}))
cWstar_tropics=dim_avg_Wrap(cWstar(:,:,{-22:22}))
eWstar_tropics=eWstar_tropics*1000.0
cWstar_tropics=cWstar_tropics*1000.0

dWstar_tropics = eWstar_tropics ; get metadata
dWstar_tropics = dWstar_tropics-cWstar_tropics


;************************************************
; create plot
;************************************************
;wks = gsn_open_wks("eps", "$IMG_PATH/WbarStarAltVsTime."+exptCase)                ; open an eps file
wks = gsn_open_wks("pdf", "$IMG_PATH/WbarStarAltVsTime."+exptCase)                ; open a pdf file
;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;gsn_define_colormap(wks,"hotres")
gsn_define_colormap(wks,"nrl_sirkes")

; Change color map using named colors
;colors = (/"white","black","lightblue","white","yellow1", \
;           "yellow3","orange1","orange2","orangered","red1","red2","red4"/)
;gsn_define_colormap(wks,colors)                  ; change colormap 

nPanel=3
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
res@cnLevelSelectionMode = "ManualLevels" 
res@cnMinLevelValF = 0.05
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 1.0
res@trYMinF = 5.
res@trYMaxF = 100.
res@gsnLeftString=" "
res@gsnRightString=" "
res@tiXAxisString       = ""
res@lbLabelBarOn         = True
res@lbOrientation        = "Vertical"
res@lbTopMarginF         =  0.0
res@lbBottomMarginF      = -3.0
res@lbRightMarginF       = -0.1
res@lbLeftMarginF        = 0.3
plot(0) = gsn_csm_pres_hgt(wks,cWstar_tropics(press|:, time|:),res)      ; contour the variable

res@lbLabelBarOn         = True
res@lbTopMarginF         = 5
res@lbBottomMarginF      = -7
res@tiXAxisString       = ""
plot(1) = gsn_csm_pres_hgt(wks,eWstar_tropics(press|:, time|:),res)      ; contour the variable

res@cnMinLevelValF = -0.5
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 0.5
res@lbLabelBarOn         = True
res@lbTopMarginF         = 5
res@lbBottomMarginF      = -7
res@tiXAxisString       = "Time (years)"
plot(2) = gsn_csm_pres_hgt(wks,dWstar_tropics(press|:, time|:),res)      ; contour the variable

resP = True
resP@gsnMaximize         = True          ; expand plot size to maximum
resP@gsnPaperOrientation = "landscape"
gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
delete(plot)  
