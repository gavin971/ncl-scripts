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
DJFpath="/h0zm.DJF/"
JJApath="/h0zm.JJA/"

testfileANN=rootpath+testcase+annpath+testcase+".cam2.h0zm.aer.annual.nc"
;testfileJJA=rootpath+testcase+JJApath+testcase+".cam2.h0zm.JJA.30-60S.nc"
;testfileDJF=rootpath+testcase+DJFpath+testcase+".cam2.h0zm.DJF.15-45N.nc"
;testsdfileJJA=rootpath+testcase+JJApath+testcase+".cam2.h0zmsd.JJA.30-60S.nc"
;testsdfileDJF=rootpath+testcase+DJFpath+testcase+".cam2.h0zmsd.DJF.15-45N.nc"
testfileJJA=rootpath+testcase+JJApath+testcase+".cam2.h0zm.JJA.60-75S.nc"
testfileDJF=rootpath+testcase+DJFpath+testcase+".cam2.h0zm.DJF.15-30N.nc"
testfileDJF2=rootpath+testcase+DJFpath+testcase+".cam2.h0zm.DJF.75-90S.nc"
testsdfileJJA=rootpath+testcase+JJApath+testcase+".cam2.h0zmsd.JJA.60-75S.nc"
testsdfileDJF=rootpath+testcase+DJFpath+testcase+".cam2.h0zmsd.DJF.15-30N.nc"
testsdfileDJF2=rootpath+testcase+DJFpath+testcase+".cam2.h0zmsd.DJF.75-90S.nc"

cntlfileANN=rootpath+cntlcase+annpath+cntlcase+".cam2.h0zm.aer.annual.nc"
;cntlfileJJA=rootpath+cntlcase+JJApath+cntlcase+".cam2.h0zm.JJA.30-60S.nc"
;cntlfileDJF=rootpath+cntlcase+DJFpath+cntlcase+".cam2.h0zm.DJF.15-45N.nc"
;cntlsdfileJJA=rootpath+cntlcase+JJApath+cntlcase+".cam2.h0zmsd.JJA.30-60S.nc"
;cntlsdfileDJF=rootpath+cntlcase+DJFpath+cntlcase+".cam2.h0zmsd.DJF.15-45N.nc"
cntlfileJJA=rootpath+cntlcase+JJApath+cntlcase+".cam2.h0zm.JJA.60-75S.nc"
cntlfileDJF=rootpath+cntlcase+DJFpath+cntlcase+".cam2.h0zm.DJF.15-30N.nc"
cntlfileDJF2=rootpath+cntlcase+DJFpath+cntlcase+".cam2.h0zm.DJF.75-90S.nc"
cntlsdfileJJA=rootpath+cntlcase+JJApath+cntlcase+".cam2.h0zmsd.JJA.60-75S.nc"
cntlsdfileDJF=rootpath+cntlcase+DJFpath+cntlcase+".cam2.h0zmsd.DJF.15-30N.nc"
cntlsdfileDJF2=rootpath+cntlcase+DJFpath+cntlcase+".cam2.h0zmsd.DJF.75-90S.nc"

testfANN=addfile(testfileANN,"r")
cntlfANN=addfile(cntlfileANN,"r")

testfJJA=addfile(testfileJJA,"r")
cntlfJJA=addfile(cntlfileJJA,"r")

testfDJF=addfile(testfileDJF,"r")
cntlfDJF=addfile(cntlfileDJF,"r")

testsdfJJA=addfile(testsdfileJJA,"r")
testsdfDJF=addfile(testsdfileDJF,"r")

cntlsdfJJA=addfile(cntlsdfileJJA,"r")
cntlsdfDJF=addfile(cntlsdfileDJF,"r")

emisfile="$MODEL_PATH/waccm/carma-soot/emissions/soot_virgin_galactic_base.nc"
emisf=addfile(emisfile,"r")

soot_emis=emisf->SOOT
soot_emis=soot_emis*1e3/365. ; Convert kg/m/y -> kg/km/day

lev_emis=emisf->lev
lev_emis=lev_emis/1e3  ; Convert m -> km

TsootJJA=testfJJA->T
print("TsootJJA:"+dimsizes(TsootJJA))
TcntlJJA=cntlfJJA->T
print("TcntlJJA:"+dimsizes(TcntlJJA))

sd1=testsdfJJA->T
print("sd1:"+dimsizes(sd1))
sd2=cntlsdfJJA->T
print("sd2:"+dimsizes(sd2))
var1=sd1(0,:)*sd1(0,:)
print("var1:"+dimsizes(var1))
var2=sd2(0,:)*sd2(0,:)
print("var2:"+dimsizes(var2))

;probtJJA = ttest(TsootJJA,var1,nyear, TcntlJJA,var2,nyear, iflag, False)
;print("probtJJA:"+dimsizes(probtJJA))

se = sqrt((var1 + var2)/sample_size)
t_JJA = abs(TsootJJA - TcntlJJA)/se

ZsootJJA=testfJJA->Z3
print("ZsootJJA:"+dimsizes(ZsootJJA))
ZcntlJJA=cntlfJJA->Z3
print("ZcntlJJA:"+dimsizes(ZcntlJJA))

TsootDJF=testfDJF->T
TcntlDJF=cntlfDJF->T
TsootsdDJF=testsdfDJF->T
print("TsootDJF:"+dimsizes(TsootDJF))
TcntlsdDJF=cntlsdfDJF->T
print("TcntlDJF:"+dimsizes(TcntlDJF))

sd1=testsdfDJF->T
print("sd1:"+dimsizes(sd1))
sd2=cntlsdfDJF->T
print("sd2:"+dimsizes(sd2))
var1=sd1(0,:)*sd1(0,:)
print("var1:"+dimsizes(var1))
var2=sd2(0,:)*sd2(0,:)
print("var2:"+dimsizes(var2))
;probtDJF = ttest(TsootDJF,var1,nyear, TcntlDJF,var2,nyear, iflag, False)
;print("probtDJF:"+dimsizes(probtDJF))
se = sqrt((var1 + var2)/sample_size)
t_DJF = abs(TsootDJF - TcntlDJF)/se

ZsootDJF=testfDJF->Z3
print("ZsootDJF:"+dimsizes(ZsootDJF))
ZcntlDJF=cntlfDJF->Z3
print("ZcntlDJF:"+dimsizes(ZcntlDJF))

; Ozone =============================
; 15-45N
O3sootJJA=testfJJA->O3
print("O3sootJJA:"+dimsizes(O3sootJJA))
O3cntlJJA=cntlfJJA->O3
print("O3cntlJJA:"+dimsizes(O3cntlJJA))

sd1=testsdfJJA->O3
print("sd1:"+dimsizes(sd1))
sd2=cntlsdfJJA->O3
print("sd2:"+dimsizes(sd2))
var1=sd1(0,:)*sd1(0,:)
print("var1:"+dimsizes(var1))
var2=sd2(0,:)*sd2(0,:)
print("var2:"+dimsizes(var2))
se = sqrt((var1 + var2)/sample_size)
t_JJA = abs(O3sootJJA - O3cntlJJA)/se

ZsootJJA=testfJJA->Z3
print("ZsootJJA:"+dimsizes(ZsootJJA))
ZcntlJJA=cntlfJJA->Z3
print("ZcntlJJA:"+dimsizes(ZcntlJJA))

; 75-90S
O3sootJJA2=testfJJA2->O3
print("O3sootJJA2:"+dimsizes(O3sootJJA2))
O3cntlJJA2=cntlfJJA2->O3
print("O3cntlJJA2:"+dimsizes(O3cntlJJA2))

sd1=testsdfJJA2->O3
print("sd1:"+dimsizes(sd1))
sd2=cntlsdfJJA2->O3
print("sd2:"+dimsizes(sd2))
var1=sd1(0,:)*sd1(0,:)
print("var1:"+dimsizes(var1))
var2=sd2(0,:)*sd2(0,:)
print("var2:"+dimsizes(var2))
se = sqrt((var1 + var2)/sample_size)
t_JJA2 = abs(O3sootJJA2 - O3cntlJJA2)/se

ZsootJJA2=testfJJA2->Z3
print("ZsootJJA2:"+dimsizes(ZsootJJA2))
ZcntlJJA2=cntlfJJA2->Z3
print("ZcntlJJA2:"+dimsizes(ZcntlJJA2))



;soot_mmr=testf->SOOT01
soot_mmr=testfANN->CARAER
soot_mmr_cntl=cntlfANN->CARAER
ZsootANN=testfANN->Z3
print("ZsootANN:"+dimsizes(ZsootANN))
ZcntlANN=cntlfANN->Z3
print("ZcntlANN:"+dimsizes(ZcntlANN))

lat=testfANN->lat
lon=testfANN->lon
lev=testfANN->lev
time=testfANN->time

nlev=dimsizes(lev)

DeltaT=new((/5,nlev/),"float")
DeltaT@_FillValue=-9999

DeltaT(0,:)=TsootDJF-TcntlDJF
DeltaT(1,:)=TsootJJA-TcntlJJA
DeltaT(2,:)=where(t_DJF.gt.sig_level,DeltaT(0,:),-9999)
DeltaT(3,:)=where(t_JJA.gt.sig_level,DeltaT(1,:),-9999)
DeltaT(4,:)=0.0

;print(TsootJJA+" | "+TcntlJJA+" | "+DeltaT(1,:))

;DeltaQRS=new((/3,nlev/),"float")
;DeltaQRS(0,:)=dim_avg_Wrap(QRSsoot(:,{45:60}))-dim_avg_Wrap(QRScntl(:,{45:60}))
;DeltaQRS(1,:)=dim_avg_Wrap(QRSsoot(:,{-45:-60}))-dim_avg_Wrap(QRScntl(:,{-45:-60}))
;;DeltaQRS(0,:)=QRSsoot(:,{45})-QRScntl(:,{45})
;;DeltaQRS(1,:)=QRSsoot(:,{-45})-QRScntl(:,{-45})
;DeltaQRS=DeltaQRS*3600.*24. ; Convert K/s -> K/day
;DeltaQRS(2,:)=0.0
;;print("lev(:)===========DeltaQRS(0,:)===========DeltaQRS(1,:)=======")
;;print (lev(:)+"   |   "+DeltaQRS(0,:)+"   |   "+DeltaQRS(1,:))

soot_mass=new((/2,nlev/),"float")
soot_mass(0,:)=dim_avg_Wrap(soot_mmr(:,{45:60}))-dim_avg_Wrap(soot_mmr_cntl(:,{45:60}))
soot_mass(1,:)=dim_avg_Wrap(soot_mmr(:,{-45:-60}))-dim_avg_Wrap(soot_mmr_cntl(:,{-45:-60}))
soot_mass=1e12*soot_mass  ; Convert kg soot/kg air -> ng soot/kg air

z2=new((/2,nlev/),"float")
;z3=new((/3,nlev/),"float")
;z4=new((/4,nlev/),"float")
z5=new((/5,nlev/),"float")

z2(0,:)=(dim_avg_Wrap(ZsootANN(:,{45:60}))+dim_avg_Wrap(ZcntlANN(:,{45:60})))/2.0
z2(1,:)=(dim_avg_Wrap(ZsootANN(:,{-45:-60}))+dim_avg_Wrap(ZcntlANN(:,{-45:-60})))/2.0
z2=z2/1e3  ; Convert m -> km

z5(0,:)=(ZsootDJF+ZcntlDJF)/2.0
z5(1,:)=(ZsootJJA+ZcntlJJA)/2.0
z5(2,:)=z5(0,:)
z5(3,:)=z5(1,:)
z5(4,:)=z5(1,:)
z5=z5/1e3  ; Convert m -> km

pdffile="$IMG_PATH/bc/bcProfilePlot.JJA60-75S.DJF15-30N"
wks   = gsn_open_wks ("pdf", pdffile)
;wks   = gsn_open_wks ("pdf", "/Volumes/Data/Images/vgProfilePlot")


  plot = new(4, graphic)

  res = True
  res@gsnDraw           = False
  res@gsnFrame          = False
  
  res@vpXF             = 0.1                  ; position and sizes
  res@vpYF             = 0.94                  ; for XY plot
  res@vpWidthF         = 0.2
  res@vpHeightF        = 0.25

  res@xyLineThicknessF = 2.5
  res@xyMonoLineThickness = True
  
  res@trYMinF = 10.0
  res@trYMaxF = 50.0
;  res@trXMinF = -0.5
  res@trXMaxF = 295.0
  res@tiYAxisString = "Altitude (km)"
  res@tiXAxisString = "BC emission rate (kg km~S~-1~N~ day~S~-1~N~)"

  plot(0) = gsn_csm_xy (wks,soot_emis,lev_emis,res)
  
  res@tiYAxisString = " "
  res@tmYLLabelsOn = False
  res@vpXF             = 0.30                  ; position and sizes
  res@trXMinF = 0.0
  res@trXMaxF = 9.0
;  delete(res@trXMinF)
;  delete(res@trXMaxF)  
  res@tiXAxisString = "BC mass (ng/kg air)"
  
  plot(1) = gsn_csm_xy (wks,soot_mass,z2,res)

  res@vpXF             = 0.50                  ; position and sizes
  res@trXMinF = -0.6
  res@trXMaxF = 1.6
  res@tiXAxisString = "~F33~D~F21~T (K)"
  res@xyMonoLineThickness = False
  res@xyLineThicknesses = (/1.5, 1.5, 5.0, 5.0, 1.0/)
  res@xyMonoDashPattern = False
  res@xyDashPatterns = (/0,1,0,1,2/)
  
  plot(2) = gsn_csm_xy (wks,DeltaT,z5,res)
  
;  res@vpXF             = 0.70                  ; position and sizes
;  res@tmYROn = True
;  res@tmYRLabelsOn = True
;  res@trXMinF = -0.018
;  res@trXMaxF = 0.018
;  res@tiXAxisString = "~F33~D~F21~Q (K day~S~-1~N~)"
  plot(4) = gsn_csm_xy (wks,DeltaO3,z3,res)

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
