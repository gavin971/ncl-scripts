load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

path="$MODEL_PATH/ccsm/run/"

test="b40.rcp4_5.2deg.wcm.carma.bc5tgIndPak.001"
cntl="b40.rcp4_5.2deg.wcm.001"

sub="h0zm"

testHistoryPattern=path+test+"/"+sub+"/"+"*.nc"
testHistoryNames = systemfunc("ls " + testHistoryPattern)
;print(testHistoryNames)
;print(dimsizes(testHistoryNames))
testHistory = addfiles(testHistoryNames, "r")
ListSetType(testHistory, "cat")

cntlHistoryPattern=path+cntl+"/"+sub+"/"+"*.nc"
cntlHistoryNames = systemfunc("ls " + cntlHistoryPattern)
;print(cntlHistoryNames)
;print(dimsizes(cntlHistoryNames))
cntlHistory = addfiles(cntlHistoryNames, "r")
ListSetType(cntlHistory, "cat")

testfile=addfile(testHistoryNames(0),"r")
gw=testfile->gw
lat=testfile->lat
nLats=dimsizes(lat)

FSNS_test  = addfiles_GetVar(testHistory, testHistoryNames, "FSNS")
FSNSC_test = addfiles_GetVar(testHistory, testHistoryNames, "FSNSC")
TS_test    = addfiles_GetVar(testHistory, testHistoryNames, "TS")
PRECT_test = addfiles_GetVar(testHistory, testHistoryNames, "PRECT")
time       = addfiles_GetVar(testHistory, testHistoryNames, "time")
time=time/365.0

FSNS_cntl  = addfiles_GetVar(cntlHistory, cntlHistoryNames, "FSNS")
FSNSC_cntl = addfiles_GetVar(cntlHistory, cntlHistoryNames, "FSNSC")
TS_cntl    = addfiles_GetVar(cntlHistory, cntlHistoryNames, "TS")
PRECT_cntl = addfiles_GetVar(cntlHistory, cntlHistoryNames, "PRECT")

TS_test_gavg = dim_avg_wgt_Wrap(TS_test,gw,0)
TS_cntl_gavg = dim_avg_wgt_Wrap(TS_cntl,gw,0)
TS_diff = TS_test_gavg - TS_cntl_gavg

;print(TS_test_gavg+" "+TS_cntl_gavg+" "+TS_diff)

PRECT_test_gavg = dim_avg_wgt_Wrap(PRECT_test,gw,0)
PRECT_cntl_gavg = dim_avg_wgt_Wrap(PRECT_cntl,gw,0)
PRECT_diff = (PRECT_test_gavg - PRECT_cntl_gavg)*1000.*3600.*24.

FSNS_test_gavg = dim_avg_wgt_Wrap(FSNS_test,gw,0)
FSNS_cntl_gavg = dim_avg_wgt_Wrap(FSNS_cntl,gw,0)
FSNS_diff = FSNS_test_gavg - FSNS_cntl_gavg

FSNSC_test_gavg = dim_avg_wgt_Wrap(FSNSC_test,gw,0)
FSNSC_cntl_gavg = dim_avg_wgt_Wrap(FSNSC_cntl,gw,0)
FSNSC_diff = FSNSC_test_gavg - FSNSC_cntl_gavg

TS_Robock = (/ (/0.072, 0.134/),\ ; (Time(y), TS change(K)) from Robock 5Tg study
               (/0.158, 0.000/),\
               (/0.230, -0.249/),\
               (/0.316, -0.537/),\
               (/0.402, -0.982/),\
               (/0.488, -1.126/),\
               (/0.646, -1.097/),\
               (/0.703, -0.963/),\
               (/0.919, -0.968/),\
               (/1.062, -0.901/),\
               (/1.148, -1.136/),\
               (/1.234, -1.284/),\
               (/1.306, -1.308/),\
               (/1.407, -1.313/),\
               (/1.507, -1.347/),\
               (/1.579, -1.313/),\
               (/1.651, -1.284/),\
               (/1.737, -1.270/),\
               (/1.852, -1.222/),\
               (/1.923, -1.141/),\
               (/1.967, -1.088/),\
               (/2.067, -1.073/),\
               (/2.167, -1.217/),\
               (/2.254, -1.284/),\
               (/2.311, -1.256/),\
               (/2.397, -1.265/),\
               (/2.498, -1.251/),\
               (/2.584, -1.256/),\
               (/2.756, -1.380/),\
               (/2.856, -1.423/),\
               (/2.914, -1.337/),\
               (/2.986, -1.112/) /)
 
PRECT_time = (/0.072,0.144,0.244,0.33,0.416,0.502,0.574,0.675,0.732,0.833,0.919,1.005,1.062,\
               1.177,1.249,1.321,1.407,1.507,1.579,1.665,1.751,1.823,1.909,1.981,2.11,2.139,\
               2.239,2.34,2.397,2.498,2.598,2.67,2.742,2.856,2.943,3.0/)
PRECT_rbck = (/-0.002,-0.019,0.006,-0.045,-0.156,-0.173,-0.235,-0.235,-0.219,-0.225,-0.257,\
               -0.204,-0.216,-0.239,-0.227,-0.276,-0.268,-0.239,-0.256,-0.245,-0.241,-0.287,\
               -0.277,-0.257,-0.245,-0.274,-0.267,-0.284,-0.292,-0.288,-0.277,-0.263,-0.248,\
               -0.29,-0.301,-0.28/)

FSNS_time = (/0.014,0.243,0.329,0.414,0.486,0.643,0.814,0.986,1.143,1.329,1.471,1.671,1.814,\
              2,2.157,2.314,2.486,2.643,2.814,2.986/)
FSNS_rbck = (/-0.339,-0.339,-6.949,-16.441,-14.746,-13.729,-13.729,-12.881,-13.39,-13.051,\
              -12.881,-12.881,-13.22,-12.542,-12.373,-11.864,-11.356,-12.034,-12.034,-11.695/)
               
  wks   = gsn_open_wks ("pdf","$IMAGE_PATH/5TgTSvsRobock") ; open workstation

  res                   = True                      ; plot mods desired
;  res@gsnPaperOrientation = "landscape"
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
;  res@vpWidthF            = 1.2           ; change aspect ratio of plot
;  res@vpHeightF           = 0.8
  res@gsnFrame     = False                        ; don't advance frame yet
  res@xyLineThicknessF = 6.0
;  res@vpWidthF=0.85

  res@xyLineColor       = "red"          ; change line color
  res@xyMonoLineColor   = True            ; same color for all lines

  res@trXMinF = 0
  res@trXMaxF = 3.0

  res@trYMaxF = 0.3
  res@trYMinF = -2.1

  res@tiYAxisString = "Surface Temperature Anomaly (~S~o~N~C)"
  res@tiXAxisString = "Year"
               
  plot = gsn_csm_xy (wks,TS_Robock(:,0),TS_Robock(:,1),res) ; draw plot
  
  res@xyLineColor       = "black"          ; change line color
  res@xyLineThicknessF = 8.0
  plot = gsn_csm_xy (wks,time,TS_diff,res) ; draw  plot
  frame(wks)
  
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.8
  res@trYMaxF = 0.0
  res@trYMinF = -0.35
  res@tiYAxisString = "Precipitation Anomaly (mm/day)"
  res@xyLineThicknessF = 6.0
  res@xyLineColor       = "red"          ; change line color
  plot = gsn_csm_xy (wks,PRECT_time,PRECT_rbck,res) ; draw plot
  res@xyLineThicknessF = 8.0
  res@xyLineColor       = "black"          ; change line color
  plot = gsn_csm_xy (wks,time,PRECT_diff,res) ; draw plot
  frame(wks)
  
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7
  res@tiYAxisString = "SW Flux Anomaly at Surface (W/m~S~2~N~)"
  res@trYMinF = -23
  res@xyLineThicknessF = 6.0
  res@xyLineColor       = "red"          ; change line color
  plot = gsn_csm_xy (wks,FSNS_time,FSNS_rbck,res) ; draw plot
  res@xyLineThicknessF = 8.0
  res@xyLineColor       = "black"          ; change line color
  plot = gsn_csm_xy (wks,time,FSNS_diff,res) ; draw plot
  res@xyLineColor       = "blue"          ; change line color
  plot = gsn_csm_xy (wks,time,FSNSC_diff,res) ; draw plot
  frame(wks)
  

