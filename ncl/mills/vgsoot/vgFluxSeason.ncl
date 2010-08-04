load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

modelpath=getenv("MODEL_PATH")

sootcase="vgsoot08"
cntlcase="vgcontrol02"

sootpath=modelpath+"/waccm/run/carma-soot/"+sootcase+"/climo/2035-2064/"
cntlpath=modelpath+"/waccm/run/carma-soot/"+cntlcase+"/climo/2035-2064/"

sootHistoryPattern=sootpath+sootcase+"_[01][0-9]_climo.nc"
sootHistoryNames = systemfunc("ls " + sootHistoryPattern)
sootHistory = addfiles(sootHistoryNames, "r")
ListSetType(sootHistory, "cat")

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

wgt = new((/ nLats, nLons /), double)
wgt = conform(wgt, gw, 0)
wgt!0="lat"
wgt&lat=lat
wgt!1="lon"
wgt&lon=lon

FLNT_Soot =addfiles_GetVar(sootHistory, sootHistoryNames, "FLNT")
FLNT_Cntl =addfiles_GetVar(cntlHistory, cntlHistoryNames, "FLNT")
FLNTC_Soot=addfiles_GetVar(sootHistory, sootHistoryNames, "FLNTC")
FLNTC_Cntl=addfiles_GetVar(cntlHistory, cntlHistoryNames, "FLNTC")

FSNT_Soot =addfiles_GetVar(sootHistory, sootHistoryNames, "FSNT")
FSNT_Cntl =addfiles_GetVar(cntlHistory, cntlHistoryNames, "FSNT")
FSNTC_Soot=addfiles_GetVar(sootHistory, sootHistoryNames, "FSNTC")
FSNTC_Cntl=addfiles_GetVar(cntlHistory, cntlHistoryNames, "FSNTC")

print("Soot")
print(dimsizes(FLNT_Soot))
;print(dimsizes(FLNTC_Soot))
;print(dimsizes(FSNT_Soot))
;print(dimsizes(FSNTC_Soot))
print("Control")
print(dimsizes(FLNT_Cntl))
;print(dimsizes(FLNTC_Cntl))
;print(dimsizes(FSNT_Cntl))
;print(dimsizes(FSNTC_Cntl))


FLNT_Soot_ave =wgt_areaave2(FLNT_Soot (time|:, lat|:, lon|:),wgt,0)
FLNTC_Soot_ave=wgt_areaave2(FLNTC_Soot(time|:, lat|:, lon|:),wgt,0)
FSNT_Soot_ave =wgt_areaave2(FSNT_Soot (time|:, lat|:, lon|:),wgt,0)
FSNTC_Soot_ave=wgt_areaave2(FSNTC_Soot(time|:, lat|:, lon|:),wgt,0)

FLNT_Cntl_ave =wgt_areaave2(FLNT_Cntl (time|:, lat|:, lon|:),wgt,0)
FLNTC_Cntl_ave=wgt_areaave2(FLNTC_Cntl(time|:, lat|:, lon|:),wgt,0)
FSNT_Cntl_ave =wgt_areaave2(FSNT_Cntl (time|:, lat|:, lon|:),wgt,0)
FSNTC_Cntl_ave=wgt_areaave2(FSNTC_Cntl(time|:, lat|:, lon|:),wgt,0)

FLNT_diff=FLNT_Soot_ave-FLNT_Cntl_ave
FLNTC_diff=FLNTC_Soot_ave-FLNTC_Cntl_ave
FSNT_diff=FSNT_Soot_ave-FSNT_Cntl_ave
FSNTC_diff=FSNTC_Soot_ave-FSNTC_Cntl_ave

print("FLNT_diff")
print(dimsizes(FLNT_diff))

do i = 0,11
  print(i+" "+FLNT_diff(i)+" "+FLNTC_diff(i)+" "+FSNT_diff(i)+" "+FSNTC_diff(i))
end do

month=ispan(0,13,1)
data=new((/4,14/),float)

data(0,1:12)=FLNT_diff
data(1,1:12)=FSNT_diff
data(2,1:12)=FLNTC_diff
data(3,1:12)=FSNTC_diff

data(:,0)=data(:,12)
data(:,13)=data(:,1)

data!1="month"
data&month=month

monthnames = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/)

;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("pdf","$IMAGE_PATH/vgFluxSeason")      ; open workstation
  res                   = True                     ; plot mods desired
  res@gsnPaperOrientation = "landscape"
  res@gsnMaximize       = True
  res@gsnBoxMargin      = 0.0
  res@xyDashPatterns = (/ 0, 0, 1, 1/)               ; Make curves all solid
  res@xyLineColors = (/"Blue","Red","Blue","Red"/)
  res@xyLineThicknessF = 3.0
  
;  res@trYLog            = True                       ; log Y-axis
  res@trYMinF           = -0.3
  res@trYMaxF           = 0.45
;  res@trYReverse        = True                       ; reverse Y-axis

;  res@trXMinF           = 0.0
;  res@trXMaxF           = 1.5
;  res@trXLog            = True                       ; log X-axis
  
 res@trXMinF = 0.5   
 res@trXMaxF = 12.5  

 res@tiYAxisString     ="Change in flux (W/m~S~2~N~)"
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
; res@gsnXYBarChart         = True            ; create bar chart 
; res@gsnAboveYRefLineColor = "blue"           ; above ref line fill red
; res@gsnBelowYRefLineColor = "red"          ; below ref line fill blue
; res@gsnAboveYRefLineBarPatterns = 1
; res@gsnBelowYRefLineBarPatterns = 1
;  res@gsnDraw               = False       ; Do not draw plot
;  res@gsnFrame              = False       ; Do not advance frame
 
;  res2=res
  
;  res@gsnAboveYRefLineColor = "blue"
;  res@gsnBelowYRefLineColor = "blue"
  
;  res2=res
 
  plot  = gsn_csm_xy (wks,month,data,res) ; create plot

