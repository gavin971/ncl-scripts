load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

testcase="bcrckt03"
cntlcase="bccntl02"

; significance level (95%)
nyear=40
sample_size = 2*nyear-2 
;sig_level = t_cvf(.025,sample_size)
sig_level = 1.99085 ; 40-year 95% significance cutoff (from IDL)
;sig_level = 1.66462
;sig_level = 1.29250
iflag= False               ; population variance similar

rootpath="$MODEL_PATH/waccm/run/carma-soot/"
annpath="/h0zm.annual/"
monpath="/h0zm.monavg/"
JJApath="/h0zm.JJA/"
DJFpath="/h0zm.DJF/"

testfileJul=rootpath+testcase+monpath+testcase+".cam2.h0zm.07.nc"
testfileJan=rootpath+testcase+monpath+testcase+".cam2.h0zm.01.nc"
testsdfileJul=rootpath+testcase+monpath+testcase+".cam2.h0zmsd.07.nc"
testsdfileJan=rootpath+testcase+monpath+testcase+".cam2.h0zmsd.01.nc"

cntlfileJul=rootpath+cntlcase+monpath+cntlcase+".cam2.h0zm.07.nc"
cntlfileJan=rootpath+cntlcase+monpath+cntlcase+".cam2.h0zm.01.nc"
cntlsdfileJul=rootpath+cntlcase+monpath+cntlcase+".cam2.h0zmsd.07.nc"
cntlsdfileJan=rootpath+cntlcase+monpath+cntlcase+".cam2.h0zmsd.01.nc"

testfileJJA=rootpath+testcase+JJApath+testcase+".cam2.h0zm.JJA.nc"
testfileDJF=rootpath+testcase+DJFpath+testcase+".cam2.h0zm.DJF.nc"
testsdfileJJA=rootpath+testcase+JJApath+testcase+".cam2.h0zm.JJA_sd.nc"
testsdfileDJF=rootpath+testcase+DJFpath+testcase+".cam2.h0zm.DJF_sd.nc"

cntlfileJJA=rootpath+cntlcase+JJApath+cntlcase+".cam2.h0zm.JJA.nc"
cntlfileDJF=rootpath+cntlcase+DJFpath+cntlcase+".cam2.h0zm.DJF.nc"
cntlsdfileJJA=rootpath+cntlcase+JJApath+cntlcase+".cam2.h0zm.JJA_sd.nc"
cntlsdfileDJF=rootpath+cntlcase+DJFpath+cntlcase+".cam2.h0zm.DJF_sd.nc"

testfJul=addfile(testfileJul,"r")
cntlfJul=addfile(cntlfileJul,"r")

testfJan=addfile(testfileJan,"r")
cntlfJan=addfile(cntlfileJan,"r")

testsdfJul=addfile(testsdfileJul,"r")
testsdfJan=addfile(testsdfileJan,"r")

cntlsdfJul=addfile(cntlsdfileJul,"r")
cntlsdfJan=addfile(cntlsdfileJan,"r")

testfJJA=addfile(testfileJJA,"r")
cntlfJJA=addfile(cntlfileJJA,"r")

testfDJF=addfile(testfileDJF,"r")
cntlfDJF=addfile(cntlfileDJF,"r")

testsdfJJA=addfile(testsdfileJJA,"r")
testsdfDJF=addfile(testsdfileDJF,"r")

cntlsdfJJA=addfile(cntlsdfileJJA,"r")
cntlsdfDJF=addfile(cntlsdfileDJF,"r")

TsootJul=testfJul->TS
print("TsootJul:"+dimsizes(TsootJul))
TcntlJul=cntlfJul->TS
print("TcntlJul:"+dimsizes(TcntlJul))
sd1=testsdfJul->TS
print("sd1:"+dimsizes(sd1))
sd2=cntlsdfJul->TS
print("sd2:"+dimsizes(sd2))
var1=sd1(0,:)*sd1(0,:)
print("var1:"+dimsizes(var1))
var2=sd2(0,:)*sd2(0,:)
print("var2:"+dimsizes(var2))
;probtJJA = ttest(TsootJJA,var1,nyear, TcntlJJA,var2,nyear, iflag, False)
;print("probtJJA:"+dimsizes(probtJJA))
se = sqrt((var1 + var2)/sample_size)
t_Jul = abs(TsootJul - TcntlJul)/se


TsootJan=testfJan->TS
TcntlJan=cntlfJan->TS
TsootsdJan=testsdfJan->TS
print("TsootJan:"+dimsizes(TsootJan))
TcntlsdJan=cntlsdfJan->TS
print("TcntlJan:"+dimsizes(TcntlJan))
sd1=testsdfJan->TS
print("sd1:"+dimsizes(sd1))
sd2=cntlsdfJan->TS
print("sd2:"+dimsizes(sd2))
var1=sd1(0,:)*sd1(0,:)
print("var1:"+dimsizes(var1))
var2=sd2(0,:)*sd2(0,:)
print("var2:"+dimsizes(var2))
;probtJan = ttest(TsootJan,var1,nyear, TcntlJan,var2,nyear, iflag, False)
;print("probtJan:"+dimsizes(probtJan))
se = sqrt((var1 + var2)/sample_size)
t_Jan = abs(TsootJan - TcntlJan)/se

lat=testfJan->lat
nlat=dimsizes(lat)

DeltaT=new((/5,nlat/),"float")
DeltaT@_FillValue=-9999

DeltaT(0,:)=TsootJan-TcntlJan
DeltaT(1,:)=TsootJul-TcntlJul
DeltaT(2,:)=where(t_Jan.gt.sig_level,DeltaT(0,:),-9999)
DeltaT(3,:)=where(t_Jul.gt.sig_level,DeltaT(1,:),-9999)
DeltaT(4,:)=0.0

;=============================================
O3colsootJJA=testfJJA->O3stcol
O3colcntlJJA=cntlfJJA->O3stcol
sd1=testsdfJJA->O3stcol
sd2=cntlsdfJJA->O3stcol
var1=sd1(0,:)*sd1(0,:)
var2=sd2(0,:)*sd2(0,:)
;probtJJA = ttest(O3stcolootJJA,var1,nyear, O3colcntlJJA,var2,nyear, iflag, False)
;print("probtJJA:"+dimsizes(probtJJA))
se = sqrt((var1 + var2)/sample_size)
t_JJA = abs(O3colsootJJA - O3colcntlJJA)/se


O3colsootDJF=testfDJF->O3stcol
O3colcntlDJF=cntlfDJF->O3stcol
O3colsootsdDJF=testsdfDJF->O3stcol
O3colcntlsdDJF=cntlsdfDJF->O3stcol
sd1=testsdfDJF->O3stcol
sd2=cntlsdfDJF->O3stcol
var1=sd1(0,:)*sd1(0,:)
var2=sd2(0,:)*sd2(0,:)
;probtDJF = ttest(O3stcolootDJF,var1,nyear, O3colcntlDJF,var2,nyear, iflag, False)
;print("probtDJF:"+dimsizes(probtDJF))
se = sqrt((var1 + var2)/sample_size)
t_DJF = abs(O3colsootDJF - O3colcntlDJF)/se

DeltaO3col=new((/5,nlat/),"float")
DeltaO3col@_FillValue=-9999

DeltaO3col(0,:)=O3colsootDJF-O3colcntlDJF
DeltaO3col(1,:)=O3colsootJJA-O3colcntlJJA
DeltaO3col(2,:)=where(t_DJF.gt.sig_level,DeltaO3col(0,:),-9999)
DeltaO3col(3,:)=where(t_JJA.gt.sig_level,DeltaO3col(1,:),-9999)
DeltaO3col(4,:)=0.0

pdffile="$IMG_PATH/bc/bcLatSeasTSO3colPlot"
wks   = gsn_open_wks ("pdf", pdffile)

  plot = new((/1,2/), graphic)

  res = True
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@gsnPaperOrientation = "portrait"
  
;  res@vpXF             = 0.01                  ; position and sizes
;  res@vpYF             = 0.94                  ; for XY plot
;  res@vpWidthF         = 0.5
  res@vpHeightF        = 0.3

  res@xyLineThicknessF = 2.5
  res@xyMonoLineThickness = True
  

  res@trYMinF = -14.0
  res@trYMaxF = 4.0
;  res@tiXAxisString = "~F33~D~F21~T (K)"
  res@xyMonoLineThickness = False
  res@xyLineThicknesses = (/1.5, 1.5, 5.0, 5.0, 1.0/)
  res@xyMonoDashPattern = False
  res@xyDashPatterns = (/0,1,0,1,2/)
  res@tmXBLabelsOn = False
  res@tmXTLabelsOn = True
  
  plot(0,0) = gsn_csm_xy (wks,lat,DeltaO3col,res)
  
  res@vpYF             = 0.50                  ; position and sizes
  res@trYMinF = -1.4
  res@trYMaxF = 2.4
;  res@tiXAxisString = "~F33~D~F21~T (K)"
  res@xyMonoLineThickness = False
  res@xyLineThicknesses = (/1.5, 1.5, 5.0, 5.0, 1.0/)
  res@xyMonoDashPattern = False
  res@xyDashPatterns = (/0,1,0,1,2/)
  res@tmXBLabelsOn = True
  res@tmXTLabelsOn = False
  
  plot(0,1) = gsn_csm_xy (wks,lat,DeltaT,res)
  
;
; Maximize plots in frame.
;
 psres = True   
; psres@gsnPaperOrientation = "landscape"                                                   
 maximize_output(wks,psres)  ; calls draw and frame for you 

;
; If you try to use gsn_panel to panel different sized plots, you
; will get unexpected results. You can turn on gsnPanelDebug to
; see what viewport values are used, and then use these values
; to set the viewport resources for each individual plot.
;  resP                      = True
;  resP@gsnMaximize          = True
;  resP@gsnPanelXWhiteSpacePercent = 0
;  resP@gsnPanelYWhiteSpacePercent = 0
;  gsn_panel(wks, plot, (/ 1, 2 /), resP)

delete(wks)

opentest = systemfunc("open "+pdffile+".pdf")
