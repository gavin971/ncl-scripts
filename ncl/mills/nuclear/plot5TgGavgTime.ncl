load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;path="$MODEL_PATH/ccsm/run/"
;path="~/csm/"
path="/data5/mmills/"

;test="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.002"
;cntl="b40.rcp4_5.2deg.wcm.carma.bc0tgIP.002"
test1="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003"
test2="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.006"
cntl1="b40.rcp4_5.2deg.wcm.carma.bc0tg.002"
cntl2="b40.rcp4_5.2deg.wcm.carma.bc0tg.003"

sub="h0gm"

print("adding files...")
testfile1=addfile(path+test1+"/"+sub+"/"+test1+".cam2.h0.globavg.nc","r")
testfile2=addfile(path+test2+"/"+sub+"/"+test2+".cam2.h0.globavg.nc","r")
cntlfile1=addfile(path+cntl1+"/"+sub+"/"+cntl1+".cam2.h0.globavg.nc","r")
cntlfile2=addfile(path+cntl2+"/"+sub+"/"+cntl2+".cam2.h0.globavg.nc","r")

gw=testfile1->gw
lat=testfile1->lat
nLats=dimsizes(lat)

print("reading FSNS_test1")
FSNS_test1  = testfile1->FSNS
print("reading FSNSC_test1")
FSNSC_test1 = testfile1->FSNSC
print("reading TS_test1")
TS_test1    = testfile1->TS
print("reading PRECT_test1")
PRECT_test1 = testfile1->PRECT
print("reading CLDTOT_test1")
CLDTOT_test1 = testfile1->CLDTOT
print("reading ASDIR_test1")
ASDIR_test1 = testfile1->ASDIR
print("reading date_x1")
date_x1     = testfile1->date
time        = testfile1->time
time=time/365.0
dx1=dimsizes(date_x1)
print("dx1:"+dx1)

print("reading FSNS_test2")
FSNS_test2  = testfile2->FSNS
print("reading FSNSC_test2")
FSNSC_test2 = testfile2->FSNSC
print("reading TS_test2")
TS_test2    = testfile2->TS
print("reading PRECT_test2")
PRECT_test2 = testfile2->PRECT
print("reading CLDTOT_test2")
CLDTOT_test2 = testfile2->CLDTOT
print("reading ASDIR_test2")
ASDIR_test2 = testfile2->ASDIR
print("reading date_x2")
date_x2     = testfile2->date
dx2=dimsizes(date_x2)
print("dx2:"+dx2)

print("reading FSNS_cntl1")
FSNS_cntl1  = cntlfile1->FSNS
print("reading FSNSC_cntl1")
FSNSC_cntl1 = cntlfile1->FSNSC
print("reading TS_cntl1")
TS_cntl1    = cntlfile1->TS
print("reading PRECT_cntl1")
PRECT_cntl1 = cntlfile1->PRECT
print("reading CLDTOT_cntl1")
CLDTOT_cntl1 = cntlfile1->CLDTOT
print("reading ASDIR_cntl1")
ASDIR_cntl1 = cntlfile1->ASDIR
print("reading date_c1")
date_c1     = cntlfile1->date
dc_1=dimsizes(date_c1)
print("dc_1:"+dc_1)

print("reading FSNS_cntl2")
FSNS_cntl2  = cntlfile2->FSNS
print("reading FSNSC_cntl2")
FSNSC_cntl2 = cntlfile2->FSNSC
print("reading TS_cntl2")
TS_cntl2    = cntlfile2->TS
print("reading PRECT_cntl2")
PRECT_cntl2 = cntlfile2->PRECT
print("reading CLDTOT_cntl2")
CLDTOT_cntl2 = cntlfile2->CLDTOT
print("reading ASDIR_cntl2")
ASDIR_cntl2 = cntlfile2->ASDIR
print("reading date_c2")
date_c2     = cntlfile2->date
dc_2=dimsizes(date_c2)
print("dc_2:"+dc_2)

dc1=ind(date_c1.eq.date_x1(0))
dc2=ind(date_c1.eq.date_x1(dx1-1))

dc3=ind(date_c2.eq.date_x2(0))
dc4=ind(date_c2.eq.date_x2(dx2-1))

print("dc1="+dc1)
print("dc2="+dc2)
print("dc3="+dc3)
print("dc4="+dc4)

print("1: "+dimsizes(TS_test1)+" "+dimsizes(TS_cntl1(dc1:dc2)))
print("2: "+dimsizes(TS_test2)+" "+dimsizes(TS_cntl2(dc3:dc4)))

print("calculating TS_diff")
TS_diff1 = TS_test1 - TS_cntl1(dc1:dc2)
TS_diff2 = TS_test2 - TS_cntl2(dc3:dc4)

print("calculating PRECT_diff")
PRECT_diff1 = (PRECT_test1 - PRECT_cntl1(dc1:dc2))*1000.*3600.*24.
PRECT_diff2 = (PRECT_test2 - PRECT_cntl2(dc3:dc4))*1000.*3600.*24.
PRECT_diff_pct1 = 100.*(PRECT_test1 - PRECT_cntl1(dc1:dc2))/PRECT_cntl1(dc1:dc2)
PRECT_diff_pct2 = 100.*(PRECT_test2 - PRECT_cntl2(dc3:dc4))/PRECT_cntl2(dc3:dc4)

print("calculating FSNS_diff")
FSNS_diff1 = FSNS_test1 - FSNS_cntl1(dc1:dc2)
;FSNS_diff2 = FSNS_test2 - FSNS_cntl2(dc3:dc4)
FSNS_diff_pct1 = 100.*(FSNS_test1 - FSNS_cntl1(dc1:dc2))/FSNS_cntl1(dc1:dc2)

print("calculating FSNSC_diff")
FSNSC_diff1 = FSNSC_test1 - FSNSC_cntl1(dc1:dc2)
;FSNSC_diff2 = FSNSC_test2 - FSNSC_cntl2(dc3:dc4)
FSNSC_diff_pct1 = 100.*(FSNSC_test1 - FSNSC_cntl1(dc1:dc2))/FSNSC_cntl1(dc1:dc2)

print("calculating CLDTOT_diff")
CLDTOT_diff1  = CLDTOT_test1 - CLDTOT_cntl1(dc1:dc2)

print("calculating ASDIR_diff")
ASDIR_diff1  = ASDIR_test1 - ASDIR_cntl1(dc1:dc2)

; (Time(y), TS change(K)) from Robock 5Tg study
TS_time = (/ 0.072,0.158,0.23,0.316,0.402,0.488,0.646,0.703,0.919,1.062,1.148,1.234,1.306,1.407,1.507,\
  1.579,1.651,1.737,1.852,1.923,1.967,2.067,2.167,2.254,2.311,2.397,2.498,2.584,2.756,2.856,2.914,2.986,\
    3.072,3.158,3.244,3.344,3.416,3.502,3.589,3.66,3.833,3.919,3.99,4.077,4.163,4.249,4.335,4.421,4.507,\
      4.579,4.665,4.751,4.823,4.909,4.995,5.081,5.182,5.254,5.34,5.426,5.512,5.584,5.67,5.742,5.842,5.928,\
      6,6.1,6.172,6.258,6.344,6.431,6.488,6.589,6.646,6.746,6.833,6.919,7.005,7.091,7.177,7.249,7.349,7.421,\
      7.507,7.593,7.679,7.766,7.837,7.923,8.01,8.096,8.182,8.268,8.354,8.426,8.512,8.612,8.684,8.756,8.842,\
      8.928,9.014,9.1,9.172,9.258,9.344,9.502,9.603,9.703,9.761,9.861,9.947,10.019,10.091,10.191,10.278,10.349 /)

TS_rbck = (/ 0.134,0,-0.249,-0.537,-0.982,-1.126,-1.097,-0.963,-0.968,-0.901,-1.136,-1.284,-1.308,\
  -1.313,-1.347,-1.313,-1.284,-1.27,-1.222,-1.141,-1.088,-1.073,-1.217,-1.284,-1.256,-1.265,-1.251,\
  -1.256,-1.38,-1.423,-1.337,-1.112,-1.097,-1.222,-1.323,-1.28,-1.16,-1.136,-1.097,-1.212,-1.26,\
  -1.188,-0.987,-0.968,-1.217,-1.179,-1.078,-0.973,-0.93,-0.887,-0.954,-1.083,-1.246,-1.155,-1.05,\
  -0.958,-0.968,-1.059,-0.925,-0.872,-0.757,-0.743,-0.863,-1.002,-1.093,-0.968,-0.676,-0.752,-0.887,\
  -0.944,-0.805,-0.666,-0.623,-0.623,-0.7,-0.8,-1.002,-0.982,-0.867,-0.824,-0.834,-0.796,-0.642,-0.565,\
  -0.494,-0.532,-0.57,-0.671,-0.772,-0.757,-0.58,-0.671,-0.815,-0.733,-0.613,-0.46,-0.417,-0.417,-0.503,\
  -0.609,-0.685,-0.724,-0.604,-0.542,-0.666,-0.618,-0.465,-0.307,-0.278,-0.398,-0.465,-0.58,-0.556,-0.503,\
  -0.527,-0.508,-0.527,-0.379 /)


PRECT_time = (/0.072,0.144,0.244,0.33,0.416,0.502,0.574,0.675,0.732,0.833,0.919,1.005,1.062,\
  1.177,1.249,1.321,1.407,1.507,1.579,1.665,1.751,1.823,1.909,1.981,2.11,2.139,2.239,2.34,\
  2.397,2.498,2.598,2.67,2.742,2.856,2.943,3,3.072,3.172,3.258,3.344,3.416,3.502,3.574,3.66,\
  3.746,3.833,3.919,3.99,4.077,4.177,4.263,4.335,4.407,4.493,4.593,4.651,4.751,4.823,4.923,\
  4.995,5.081,5.167,5.254,5.34,5.426,5.512,5.598,5.67,5.742,5.842,5.928,6,6.086,6.172,6.258,\
  6.344,6.431,6.502,6.589,6.675,6.761,6.833,6.933,7.005,7.091,7.163,7.249,7.335,7.421,7.493,\
  7.593,7.679,7.751,7.837,7.923,8.01,8.096,8.182,8.254,8.34,8.426,8.512,8.584,8.67,8.756,8.856,\
  8.943,9.014,9.086,9.172,9.258,9.344,9.416,9.545,9.603,9.703,9.761,9.861,9.947,10.005,10.105,\
  10.177,10.263,10.364/)
  
PRECT_rbck = (/-0.002,-0.019,0.006,-0.045,-0.156,-0.173,-0.235,-0.235,-0.219,-0.225,-0.257,-0.204,-0.216,\
  -0.239,-0.227,-0.276,-0.268,-0.239,-0.256,-0.245,-0.241,-0.287,-0.277,-0.257,-0.245,-0.274,-0.267,\
  -0.284,-0.292,-0.288,-0.277,-0.263,-0.248,-0.29,-0.301,-0.28,-0.262,-0.277,-0.271,-0.283,-0.286,\
  -0.292,-0.269,-0.249,-0.259,-0.265,-0.308,-0.264,-0.261,-0.264,-0.25,-0.256,-0.239,-0.233,-0.251,\
  -0.215,-0.225,-0.245,-0.249,-0.23,-0.21,-0.216,-0.205,-0.224,-0.201,-0.206,-0.213,-0.197,-0.194,\
  -0.211,-0.198,-0.184,-0.162,-0.156,-0.156,-0.177,-0.158,-0.166,-0.177,-0.138,-0.145,-0.158,-0.174,\
  -0.128,-0.111,-0.12,-0.124,-0.127,-0.115,-0.123,-0.125,-0.106,-0.099,-0.121,-0.119,-0.075,-0.094,\
  -0.099,-0.073,-0.103,-0.097,-0.076,-0.089,-0.073,-0.07,-0.072,-0.095,-0.044,-0.071,-0.05,-0.054,\
  -0.074,-0.068,-0.064,-0.077,-0.075,-0.031,-0.07,-0.054,-0.031,-0.03,-0.037,-0.041,-0.031/)
  

FSNS_time = (/0.014,0.243,0.329,0.414,0.486,0.643,0.814,0.986,1.143,1.329,1.471,1.671,1.814,2,2.157,2.314,\
  2.486,2.643,2.814,2.986,3.129,3.3,3.471,3.643,3.8,3.986,4.143,4.3,4.471,4.657,4.814,4.986,5.143,5.3,\
  5.471,5.629,5.8,5.971,6.129,6.3,6.486,6.643,6.786,6.971,7.143,7.3,7.529,7.757,7.957,8.143,8.314,\
  8.457,8.629,8.786,8.971,9.129,9.3,9.471,9.629,9.8,9.971,10.129,10.286/)
FSNS_rbck = (/-0.339,-0.339,-6.949,-16.441,-14.746,-13.729,-13.729,-12.881,-13.39,-13.051,-12.881,-12.881,-13.22,\
  -12.542,-12.373,-11.864,-11.356,-12.034,-12.034,-11.695,-10.847,-10.339,-10,-10.339,-10.339,-10,-9.492,-8.814,\
  -8.305,-8.644,-8.644,-8.136,-7.797,-6.949,-6.78,-6.78,-6.78,-6.441,-5.932,-5.424,-5.254,-5.254,-5.254,-5.085,\
  -4.576,-4.068,-3.559,-3.559,-3.39,-3.051,-3.051,-2.712,-2.712,-2.712,-2.373,-2.203,-2.034,-2.034,-1.864,-1.864,\
  -1.525,-1.525,-1.186/)
               
  wks   = gsn_open_wks ("pdf","$IMAGE_PATH/TSvsRobock."+test1+"."+test2) ; open workstation

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

  res@tmXBMaxTicks     = 15
  res@tmXBMinorPerMajor= 11
  res@tmXBAutoPrecision= False
  res@tmXBPrecision    = 2
  
  res@trXMinF = 0
  res@trXMaxF = max(time)

  res@trYMaxF = 0.2
  res@trYMinF = -2.8

  res@tiYAxisString = "Surface Temperature Anomaly (~S~o~N~C)"
  res@tiXAxisString = "Year"
               
  plot = gsn_csm_xy (wks,TS_time,TS_rbck,res) ; draw plot
  
  res@xyLineColor       = "black"          ; change line color
  res@xyLineThicknessF = 8.0
  plot = gsn_csm_xy (wks,time,TS_diff1,res) ; draw  plot
  res@xyLineColor       = "green"          ; change line color
  res@xyLineThicknessF = 8.0
  plot = gsn_csm_xy (wks,time,TS_diff2,res) ; draw  plot
  frame(wks)
  
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.8
  res@trYMaxF = 0.0
  res@trYMinF = -0.4
  res@tiYAxisString = "Precipitation Anomaly (mm/day)"
  res@xyLineThicknessF = 6.0
  res@xyLineColor       = "red"          ; change line color
  plot = gsn_csm_xy (wks,PRECT_time,PRECT_rbck,res) ; draw plot
  res@xyLineThicknessF = 8.0
  res@xyLineColor       = "black"          ; change line color
  plot = gsn_csm_xy (wks,time,PRECT_diff1,res) ; draw plot
  res@xyLineThicknessF = 8.0
  res@xyLineColor       = "green"          ; change line color
  plot = gsn_csm_xy (wks,time,PRECT_diff2,res) ; draw plot
  frame(wks)
  
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.8
  res@trYMaxF = 0.0
  res@trYMinF = -15
  res@tiYAxisString = "Precipitation Anomaly (%)"
  res@xyLineThicknessF = 8.0
  res@xyLineColor       = "black"          ; change line color
  plot = gsn_csm_xy (wks,time,PRECT_diff_pct1,res) ; draw plot
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
  plot = gsn_csm_xy (wks,time,FSNS_diff1,res) ; draw plot
  res@xyLineThicknessF = 8.0
;  res@xyLineColor       = "green"          ; change line color
;  plot = gsn_csm_xy (wks,time,FSNS_diff2,res) ; draw plot
  res@xyLineColor       = "blue"          ; change line color
  plot = gsn_csm_xy (wks,time,FSNSC_diff1,res) ; draw plot
;  res@xyLineColor       = "maroon"          ; change line color
;  plot = gsn_csm_xy (wks,time,FSNSC_diff2,res) ; draw plot
  frame(wks)
  
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7
  res@tiYAxisString = "SW Flux Anomaly at Surface (%)"
  res@trYMinF = -15
  res@xyLineThicknessF = 8.0
  res@xyLineColor       = "black"          ; change line color
  plot = gsn_csm_xy (wks,time,FSNS_diff_pct1,res) ; draw plot
  res@xyLineColor       = "blue"          ; change line color
  plot = gsn_csm_xy (wks,time,FSNSC_diff_pct1,res) ; draw plot
  frame(wks)
  
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7
  res@tiYAxisString = "CLDTOT (frac)"
  res@trYMinF=0.52
  res@trYMaxF=0.6
  res@xyLineThicknessF = 8.0
  res@xyLineColor       = "black"          ; change line color
  plot = gsn_csm_xy (wks,time,CLDTOT_test1,res) ; draw plot
  res@xyLineColor       = "blue"          ; change line color
  plot = gsn_csm_xy (wks,time,CLDTOT_cntl1,res) ; draw plot
;  res@xyLineColor       = "red"          ; change line color
;  plot = gsn_csm_xy (wks,time,CLDTOT_diff1,res) ; draw plot
  frame(wks)
  
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7
  res@tiYAxisString = "ASDIR (frac)"
  res@trYMinF=0.55
  res@trYMaxF=0.59
  res@xyLineThicknessF = 8.0
  res@xyLineColor       = "black"          ; change line color
  plot = gsn_csm_xy (wks,time,ASDIR_test1,res) ; draw plot
  res@xyLineColor       = "blue"          ; change line color
  plot = gsn_csm_xy (wks,time,ASDIR_cntl1,res) ; draw plot
;  res@xyLineColor       = "red"          ; change line color
;  plot = gsn_csm_xy (wks,time,ASDIR_diff1,res) ; draw plot
  frame(wks)
  

