;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "../common.ncl"

path="$MODEL_PATH/ccsm/case/"

exptCase="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"
cntlCase="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"

exptFilepath=path+exptCase+"/proc/concat_tem_diag0_"+exptCase+".nc"
cntlFilepath=path+cntlCase+"/proc/concat_tem_diag0_"+cntlCase+".nc"

exptFile=addfile(exptFilepath,"r")
cntlFile=addfile(cntlFilepath,"r")

eWstar = exptFile->WSTAR
cWstar = cntlFile->WSTAR

press = exptFile->lev
press@units="hPa"
press@long_name="Pressure"
nlev=dimsizes(press)

;eWstar!1="press"
;eWstar&press=press
;eWstar&time=time

;cWstar!1="press"
;cWstar&press=press
;cWstar&time=time

eWstar_tropics=dim_avg_Wrap(eWstar(:,:,{-22:22}))
cWstar_tropics=dim_avg_Wrap(cWstar(:,:,{-22:22}))
eWstar_tropics=eWstar_tropics*1000.0
cWstar_tropics=cWstar_tropics*1000.0

;ntimes=9*12
ntimes=10
eWstar_ann=new((/ntimes,nlev/),typeof(eWstar))
eWstar_ann!1="press"
eWstar_ann&press=press
eWstar_ann!0="time"
time=fspan(0.5,9.5,ntimes)
time@units="years"
;printVarSummary(eWstar_ann)
;printVarSummary(time)
eWstar_ann&time=time

cWstar_ann=eWstar_ann
do i = 0,ntimes-1
  if (ntimes.eq.10) then
    i1=i*12
    i2=(i+1)*12-1
  else ; boxcar smooth
    i1=i
    i2=i+11
  end if
;  print(i+" "+i2+" "+9*12)
  cWstar_ann(i,:)=dim_avg_n(cWstar_tropics(i1:i2,:),0)
  eWstar_ann(i,:)=dim_avg_n(eWstar_tropics(i1:i2,:),0)
end do


dWstar_ann = eWstar_ann ; get metadata
dWstar_ann = dWstar_ann-cWstar_ann


;************************************************
; create plot
;************************************************
;wks = gsn_open_wks("eps", "$IMG_PATH/WbarStarAltVsTime."+exptCase)                ; open an eps file
wks = gsn_open_wks("pdf", "$IMG_PATH/WbarStarAltVsTime."+exptCase)                ; open a pdf file
;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
;gsn_define_colormap(wks,"hotres")
gsn_define_colormap(wks,"nrl_sirkes")
;gsn_define_colormap(wks,"WhViBlGrYeOrReWh")

; Change color map using named colors
;colors = (/"white","black","lightblue","white","yellow1", \
;           "yellow3","orange1","orange2","orangered","red1","red2","red4"/)
;gsn_define_colormap(wks,colors)                  ; change colormap 

nPanel=3
plot = new(nPanel, graphic)

res                       = True     ; plot mods desired
;res@vpHeightF             = 0.82 / nPanel    ; change aspect ratio
res@vpHeightF             = 0.62 / nPanel    ; change aspect ratio
;res@vpWidthF              = 1.0              ; change aspect ratio
res@gsnDraw = False
res@gsnFrame = False
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = True     ; turn on contour lines
res@cnLineLabelsOn        = True     ; turn on contour line labels
res@cnLineLabelFontHeightF=0.01
res@cnLineLabelFormat="@f"
res@tmYLMode    = "Explicit" 
res@tmYLValues  = (/    1E-5 ,   1E-4  , 0.001 , 0.01 , 0.1 , 1., 10., 100./)
res@tmYLLabels  = (/"10~S~-5","10~S~-4","0.001","0.01","0.1","1","10","100"/)
res@tmXBMode    = "Manual" 
res@tmXBTickSpacingF  = 1.0
res@tmXBTickStartF = 1.0
res@cnLevelSelectionMode = "ManualLevels" 
res@cnMinLevelValF = 0.05
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 0.85
res@trYMinF = 4.
res@trYMaxF = 110.
res@gsnLeftString=" "
res@gsnRightString=" "
res@tiXAxisString       = ""
res@lbLabelBarOn         = True
res@lbOrientation        = "Vertical"
res@lbTopMarginF         =  0.0
res@lbBottomMarginF      = -1.6
res@lbRightMarginF       = 0.03
res@lbLeftMarginF        = 0.30
res@lbLabelFontHeightF   = 0.01
res@tiYAxisString        = "Pressure (mb)"
res@tiYAxisFontHeightF   = 0.014
res@tiXAxisFontHeightF   = 0.014
;res@lbTitleOn            = True
;res@lbTitleString        = "mm/s"
;res@lbTitlePosition      = "Right"
;res@tmXBMajorLengthF     = 0.005
;res@tmXBMinorLengthF     = 0.0
res@tmXBMinorOn          = False
res@tmYLMinorOn          = True
res@tmXBAutoPrecision    = False
res@tmXBPrecision        = 1
res@vpYF                 = 0.95
;res@tmXBMajorLengthF     = 0.0
;res@tmXTMajorLengthF     = 0.05
;res@tmYLMajorLengthF     = 0.007
;res@tmYRMajorLengthF     = 0.007
res@tmXBOn       = True
res@tmXBLabelsOn = False
res@tmXTOn       = True
res@tmXTLabelsOn = True
res@tiYAxisOn    = False
plot(0) = gsn_csm_pres_hgt(wks,cWstar_ann(press|:, time|:),res)      ; contour the variable

;==================================================================================
res@lbLabelBarOn         = True
res@lbTopMarginF         = 7
res@lbBottomMarginF      = -7
res@tiXAxisString        = ""
;res@tmXBMajorLengthF     = 0.0
;res@tmXTMajorLengthF     = 0.0
;res@vpYF                 = 0.675
res@vpYF                 = 0.743
res@tmXBOn       = True
res@tmXBLabelsOn = False
res@tmXTOn       = True
res@tmXTLabelsOn = False
res@tiYAxisOn    = True
plot(1) = gsn_csm_pres_hgt(wks,eWstar_ann(press|:, time|:),res)      ; contour the variable

;==================================================================================
;res@cnMinLevelValF = -0.4
;res@cnLevelSpacingF = 0.05
;res@cnMaxLevelValF = 0.4
res@cnLevelSelectionMode = "ExplicitLevels" 
delete(res@cnLevelSpacingF)
res@cnLevels =(/-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,.15,.2,.25,.3,.35/)
res@lbLabelBarOn         = True
res@lbRightMarginF       = -0.03
res@lbLeftMarginF        = 0.33
res@lbTopMarginF         = 0.0
res@lbBottomMarginF      = -0.3
res@tiXAxisString        = "Time (years)"
;delete(res@tmXBMajorLengthF)
;res@vpYF                 = 0.4
res@vpYF                 = 0.535
res@tmXBOn       = True
res@tmXBLabelsOn = True
res@tmXTOn       = True
res@tmXTLabelsOn = False
res@tiYAxisOn    = False
plot(2) = gsn_csm_pres_hgt(wks,dWstar_ann(press|:, time|:),res)      ; contour the variable

;==================================================================================
;frame(wks)
;resP = True
;resP@gsnMaximize         = True          ; expand plot size to maximum
;resP@gsnPaperOrientation = "landscape"
;gsn_panel(wks, plot, (/ nPanel, 1 /), resP)
; See what happens in PS file when you maximize all the plots.
psres = True   
maximize_output(wks,psres) 
;draw(wks)
