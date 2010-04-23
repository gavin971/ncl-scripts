load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;modelpath=getenv("MODEL_PATH")

sootcase="vgsoot08"
cntlcase="vgcontrol02"

;runpath=modelpath+"/waccm/run/carma-soot/"
runpath="$MODEL_PATH/waccm/run/carma-soot/"

sootpath=runpath+sootcase+"/climo/2035-2064/"
sootHistoryPattern=sootpath+sootcase+"_[01][0-9]_climo.nc"
sootHistoryNames = systemfunc("ls " + sootHistoryPattern)
sootHistory = addfiles(sootHistoryNames, "r")
ListSetType(sootHistory, "cat")

cntlpath=runpath+cntlcase+"/climo/2035-2064/"
cntlHistoryPattern=cntlpath+cntlcase+"_[01][0-9]_climo.nc"
cntlHistoryNames = systemfunc("ls " + cntlHistoryPattern)
cntlHistory = addfiles(cntlHistoryNames, "r")
ListSetType(cntlHistory, "cat")

sootfile=addfile(sootHistoryNames(0),"r")
gw=sootfile->gw
lat=sootfile->lat
lon=sootfile->lon
nLats=dimsizes(lat)
nLons=dimsizes(lon)

icefracSoot=addfiles_GetVar(sootHistory, sootHistoryNames, "ICEFRAC")
icefracCntl=addfiles_GetVar(cntlHistory, cntlHistoryNames, "ICEFRAC")

print("Soot")
print(dimsizes(icefracSoot))
;print("Soot lat")
;print(icefracSoot&lat)
print("Control")
print(dimsizes(icefracCntl))
;print("Control lat")
;print(icefracCntl&lat)

wgt = new((/ nLats, nLons /), double)
wgt = conform(wgt, gw, 0)
wgt!0="lat"
wgt&lat=lat
wgt!1="lon"
wgt&lon=lon

print("Soot")
print(dimsizes(icefracSoot(:,{1:90},:)))
print("wgt")
print(dimsizes(wgt({1:90},:)))

icefracSootNH = wgt_areasum2(icefracSoot(:,{1:91},:), wgt({1:91},:), 0)
icefracSootSH = wgt_areasum2(icefracSoot(:,{-1:-91},:), wgt({-1:-91},:), 0)

print("Control")
print(dimsizes(icefracCntl(:,{1:90},:)))

icefracCntlNH = wgt_areasum2(icefracCntl(:,{1:91},:), wgt({1:91},:), 0)
icefracCntlSH = wgt_areasum2(icefracCntl(:,{-1:-91},:), wgt({-1:-91},:), 0)

pctNH=100.0*(icefracSootNH-icefracCntlNH)/icefracCntlNH
pctSH=100.0*(icefracSootSH-icefracCntlSH)/icefracCntlSH

print("NH")
do i = 0,11
  m=i+1
  print(m+" "+icefracSootNH(i)+" "+icefracCntlNH(i)+" "+pctNH(i))
end do

icefracCntlNHavg=dim_avg(icefracCntlNH)
icefracSootNHavg=dim_avg(icefracSootNH)
pctNHave=(icefracSootNHavg-icefracCntlNHavg)/icefracCntlNHavg
print("average:"+pctNHave)

print("SH")
do i = 0,11
  m=i+1
  print(m+" "+icefracSootSH(i)+" "+icefracCntlSH(i)+" "+pctSH(i))
end do

icefracCntlSHavg=dim_avg(icefracCntlSH)
icefracSootSHavg=dim_avg(icefracSootSH)
pctSHave=(icefracSootSHavg-icefracCntlSHavg)/icefracCntlSHavg
print("average:"+pctSHave)

month=ispan(0,13,1)
data=new(14,float)
data2=new(14,float)
data2(1:12)=pctNH
data(1:12)=pctSH

data(0)=data(12)
data(13)=data(1)

data2(0)=data(12)
data2(13)=data(1)

data!0="month"
data&month=month

data2!0="month"
data2&month=month

monthnames = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/)

;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("pdf","$IMAGE_PATH/vgIceFracChange")      ; open workstation
  res                   = True                     ; plot mods desired
  res@gsnPaperOrientation = "landscape"
  res@gsnMaximize       = True
  res@gsnBoxMargin      = 0.0
;  res@xyDashPattern = 0               ; Make curves all solid
;  res@xyLineColors = (/"Red","Blue"/)
  
;  res@trYLog            = True                       ; log Y-axis
  res@trYMinF           = -20
  res@trYMaxF           = 2.0
;  res@trYReverse        = True                       ; reverse Y-axis

;  res@trXMinF           = 0.0
;  res@trXMaxF           = 1.5
;  res@trXLog            = True                       ; log X-axis
  
 res@trXMinF = 0.5   
 res@trXMaxF = 12.5  

 res@tiYAxisString     ="Change in hemispheric sea ice extent (%)"
; res@tiXAxisString     ="Month"
;  res@tiMainString      = "Tropics Average"             ; add title
 res@tiYAxisFontHeightF = 0.018

 res@tmXBLabels = monthnames
 res@tmXBMode           = "Explicit"         ; Define own tick mark labels.
 res@tmXBValues         = data&month(1:12)
 res@tmXBMinorValues    = data&month(0:12)+0.5
 res@tmXBLabelFontHeightF = 0.015
 res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
 res@tmXBMinorLengthF        = 0.02    ; minor tickmarks
 res@tmXBMinorThicknessF     = 2.0
 res@tmXBMajorLengthF        = 0.00
 res@tmYLMajorOutwardLengthF        = 0.02
 res@tmYLMajorLengthF        = 0.02
  
 res@gsnYRefLine           = 0.              ; reference line   
 res@gsnXYBarChart         = True            ; create bar chart 
; res@gsnAboveYRefLineColor = "blue"           ; above ref line fill red
; res@gsnBelowYRefLineColor = "red"          ; below ref line fill blue
; res@gsnAboveYRefLineBarPatterns = 1
; res@gsnBelowYRefLineBarPatterns = 1
  res@gsnDraw               = False       ; Do not draw plot
  res@gsnFrame              = False       ; Do not advance frame
 
;  res2=res
  
  res@gsnAboveYRefLineColor = "blue"
  res@gsnBelowYRefLineColor = "blue"
  
  res2=res
 
  plot  = gsn_csm_xy (wks,month,data,res) ; create plot
  
  res2@gsnAboveYRefLineColor = "red"
  res2@gsnBelowYRefLineColor = "red"
;  res2@gsnAboveYRefLineBarPatterns = 14
;  res2@gsnBelowYRefLineBarPatterns = 14
  plot2  = gsn_csm_xy (wks,month,data2,res2) ; create plot
  
  overlay (plot, plot2)
  draw(plot)
  frame(wks)
