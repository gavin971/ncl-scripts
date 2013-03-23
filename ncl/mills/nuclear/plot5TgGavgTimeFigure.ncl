load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

host   = systemfunc("hostname")
ostype = systemfunc("echo $OSTYPE")
MODEL_PATH=getenv("MODEL_PATH")
IMAGE_PATH=getenv("IMG_PATH")

if (host.eq."modeling1") then
  path="/data4/mmills/"
else if (isStrSubset(host,"caldera").or.isStrSubset(host,"geyser").or.isStrSubset(host,"yslogin")) then
  path=MODEL_PATH+"/"
else
  path=MODEL_PATH+"/cesm/case/"  
end if
end if
;path="~/csm/"
print("path="+path)

root="b40.rcp4_5.2deg.wcm.carma."
sub="h0gm"

test1="bc5tgIP"
cntl1="bc0tg"

massfilepath=path+root+test1+".avg/proc/"

massfile=root+test1+".mass.avg.nc"
print("reading "+massfilepath+massfile)
mfileid = addfile(massfilepath+massfile, "r")
dmass = mfileid->SOOT01
delete(mfileid)
nTimes = dimsizes(dmass)

dmass_range = new((/2,nTimes/), typeof(dmass))

massfile=root+test1+".mass.max.nc"
print("reading "+massfilepath+massfile)
mfileid = addfile(massfilepath+massfile, "r")
dmass_range(0,0:nTimes-1) = mfileid->SOOT01(0:nTimes-1)
delete(mfileid)

massfile=root+test1+".mass.min.nc"
print("reading "+massfilepath+massfile)
mfileid = addfile(massfilepath+massfile, "r")
dmass_range(1,0:nTimes-1) = mfileid->SOOT01(0:nTimes-1)
delete(mfileid)
      
avgfile=path+root+test1+".avg/"+sub+"/"+root+test1+"-"+cntl1+".cam.h0.avg.nc"
maxfile=path+root+test1+".avg/"+sub+"/"+root+test1+"-"+cntl1+".cam.h0.max.nc"
minfile=path+root+test1+".avg/"+sub+"/"+root+test1+"-"+cntl1+".cam.h0.min.nc"
print("adding avg file: "+avgfile)
avgfileh=addfile(avgfile,"r")
print("adding max file: "+maxfile)
maxfileh=addfile(maxfile,"r")
print("adding min file: "+minfile)
minfileh=addfile(minfile,"r")

gw   = avgfileh->gw
lat  = avgfileh->lat
time = avgfileh->time

time   = time/365.0
nLats  = dimsizes(lat)
nTimes = dimsizes(time)
print("nTimes="+nTimes)

; Net shortwave flux at the surface
print("reading FSNS_diff")
FSNS_diff = avgfileh->FSNS(0:nTimes-1)
FSNS_range = new((/2,nTimes/), typeof(FSNS_diff))
print("reading FSNS_max")
FSNS_range(0,0:nTimes-1) = maxfileh->FSNS(0:nTimes-1)
print("reading FSNS_min")
FSNS_range(1,0:nTimes-1) = minfileh->FSNS(0:nTimes-1)

; Surface temperature
print("reading TS_diff")
TS_diff = avgfileh->TS(0:nTimes-1)
TS_range = new((/2,nTimes/), typeof(TS_diff))
print("reading TS_max")
TS_range(0,0:nTimes-1) = maxfileh->TS(0:nTimes-1)
print("reading TS_min")
TS_range(1,0:nTimes-1) = minfileh->TS(0:nTimes-1)

; Total precipitation
print("reading PRECT_diff")
PRECT_diff = avgfileh->PRECT(0:nTimes-1)
PRECT_range = new((/2,nTimes/), typeof(PRECT_diff))
print("reading PRECT_max")
PRECT_range(0,0:nTimes-1) = maxfileh->PRECT(0:nTimes-1)
print("reading PRECT_min")
PRECT_range(1,0:nTimes-1) = minfileh->PRECT(0:nTimes-1)

; convert m/s to mm/day
PRECT_diff = PRECT_diff * 24. * 3600. * 1000.
PRECT_range = PRECT_range * 24. * 3600. * 1000.

; Global mass burden from Robock 5Tg study
robock = (/ (/0.353, 4.948/), (/0.409, 4.892/), (/0.465, 4.837/), (/0.539, 4.795/), (/0.613, 4.768/), \
            (/0.687, 4.754/), (/0.761, 4.740/), (/0.873, 4.698/), (/0.984, 4.657/), (/1.188, 4.587/), \
            (/1.410, 4.504/), (/1.707, 4.421/), (/2.115, 4.310/), (/2.486, 4.129/), (/2.764, 4.005/), \
            (/3.060, 3.880/), (/3.468, 3.616/), (/3.783, 3.422/), (/4.099, 3.283/), (/4.358, 3.075/), \
            (/4.803, 2.798/), (/5.063, 2.673/), (/5.359, 2.493/), (/5.749, 2.243/), (/6.027, 2.118/), \
            (/6.305, 1.938/), (/6.583, 1.771/), (/6.787, 1.660/), (/7.139, 1.535/), (/7.510, 1.341/), \
            (/7.862, 1.202/), (/8.103, 1.119/), (/8.363, 0.994/), (/8.789, 0.856/), (/9.067, 0.772/), \
            (/9.364, 0.689/), (/9.568, 0.634/), (/9.791, 0.578/), (/9.939, 0.550/), (/10.050, 0.537/) /)

robock(:,0)=robock(:,0)-134./365. ; Robock starts May 15, vs Jan 1

; (Time(y), TS change(K)) from Robock 5Tg study
TS_time = (/ 0.072,0.158,0.23,0.316,0.402,0.488,0.646,0.703,0.919,1.062,1.148,1.234,1.306,1.407,1.507,\
  1.579,1.651,1.737,1.852,1.923,1.967,2.067,2.167,2.254,2.311,2.397,2.498,2.584,2.756,2.856,2.914,2.986,\
    3.072,3.158,3.244,3.344,3.416,3.502,3.589,3.66,3.833,3.919,3.99,4.077,4.163,4.249,4.335,4.421,4.507,\
      4.579,4.665,4.751,4.823,4.909,4.995,5.081,5.182,5.254,5.34,5.426,5.512,5.584,5.67,5.742,5.842,5.928,\
      6,6.1,6.172,6.258,6.344,6.431,6.488,6.589,6.646,6.746,6.833,6.919,7.005,7.091,7.177,7.249,7.349,7.421,\
      7.507,7.593,7.679,7.766,7.837,7.923,8.01,8.096,8.182,8.268,8.354,8.426,8.512,8.612,8.684,8.756,8.842,\
      8.928,9.014,9.1,9.172,9.258,9.344,9.502,9.603,9.703,9.761,9.861,9.947,10.019,10.091,10.191,10.278,10.349 /)

TS_time=TS_time-134./365. ; Robock starts May 15, vs Jan 1

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
  
PRECT_time=PRECT_time-134./365. ; Robock starts May 15, vs Jan 1
  
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
  
FSNS_time=FSNS_time-134./365. ; Robock starts May 15, vs Jan 1
  
FSNS_rbck = (/-0.339,-0.339,-6.949,-16.441,-14.746,-13.729,-13.729,-12.881,-13.39,-13.051,-12.881,-12.881,-13.22,\
  -12.542,-12.373,-11.864,-11.356,-12.034,-12.034,-11.695,-10.847,-10.339,-10,-10.339,-10.339,-10,-9.492,-8.814,\
  -8.305,-8.644,-8.644,-8.136,-7.797,-6.949,-6.78,-6.78,-6.78,-6.441,-5.932,-5.424,-5.254,-5.254,-5.254,-5.085,\
  -4.576,-4.068,-3.559,-3.559,-3.39,-3.051,-3.051,-2.712,-2.712,-2.712,-2.373,-2.203,-2.034,-2.034,-1.864,-1.864,\
  -1.525,-1.525,-1.186/)
 
  pdffile="TSvsRobockFigure."+test1+"-"+cntl1
  wks   = gsn_open_wks ("pdf",IMAGE_PATH+"/"+pdffile) ; open workstation

  nl = integertochar(10)  ; newline character

  ; Create a panel plot
  resP = True
  nPanel=3
  plot = new(nPanel, graphic)
  
  resP@txString = "Years after conflict"
  resP@gsnAttachPlotsXAxis =True
  resP@gsnPaperOrientation = "landscape"
  resP@gsnMaximize         = True          ; expand plot size to maximum
 
  
  res                   = True                      ; plot mods desired
  res@vpWidthF            = 1.3           ; change aspect ratio of plot
  res@vpHeightF           = 0.3
  res@gsnDraw = False
  res@gsnFrame     = False                        ; don't advance frame yet
  res@xyLineThicknessF = 5.0
;  res@vpWidthF=0.85

  res@xyLineColor       = "red"          ; change line color
  res@xyMonoLineColor   = True            ; same color for all lines

  res@tmXBMaxTicks     = 15
  res@tmXBMinorPerMajor= 3
  res@tmXTMinorPerMajor= 3
  res@tmXBAutoPrecision= False
  res@tmXBPrecision    = 2
  
  res@trXMinF = 0.0
  res@trXMaxF = max(time)
  res@trYMinF = 0.0
  res@trYMaxF = 5.0
  res@tiXAxisString = "Year"
;  res@tmXTLabelsOn = True
;  res@tmYRLabelsOn = True

  res@trYMinF = -0.01
;  res@trYMaxF = 5.0
  res@tmYLPrecision = 1
  res@tmYLAutoPrecision = False
  res@tiYAxisString = "Mass burden (Tg)"
  
  res2=res
  res2@xyLineThicknessF = 0.5
  res2@xyMonoDashPattern = True
  res2@xyLineColor       = "lightblue"
  res2@gsnXYFillColors = (/"lightblue"/)

  res3=res
  res3@xyLineThicknessF = 4.0
  res3@xyLineColor       = "midnightblue"
  
  baseplot = gsn_csm_xy (wks, dmass&time/365., dmass_range/1e9,res2) ; draw plot
  overlay(baseplot, gsn_csm_xy (wks, robock(:,0), robock(:,1), res)) ; draw  plot
  overlay(baseplot, gsn_csm_xy (wks, dmass&time/365., dmass/1e9,res3)) ; draw  plot
  
  res@tiYAxisString = "SW Flux (W m~S~-2~N~)"
;  res@trYMaxF = -0.01
  res@trYMaxF = 2.0
  res@trYMinF = -17.0
  res@tmYLPrecision=2
  
  res2=res
  res2@xyLineThicknessF = 0.5
  res2@xyMonoDashPattern = True
  res2@xyLineColor       = "lightblue"
  res2@gsnXYFillColors = (/"lightblue"/)
  res3=res
  res3@xyLineThicknessF = 3.0
  res3@xyLineColor       = "midnightblue"
  
  plot(0) = gsn_csm_xy (wks,time,FSNS_range,res2) ; draw plot
  overlay(plot(0), gsn_csm_xy (wks,FSNS_time,FSNS_rbck,res)) ; draw plot
  overlay(plot(0), gsn_csm_xy (wks,time,FSNS_diff,res3)) ; draw plot
  
  res@tiYAxisString = "Temperature (K)"
  res@trYMaxF = 0.1
  res@trYMinF = -2.1

  res2=res
  res2@xyLineThicknessF = 0.5
  res2@xyMonoDashPattern = True
  res2@xyLineColor       = "lightblue"
  res2@gsnXYFillColors = (/"lightblue"/)
  res3=res
  res3@xyLineThicknessF = 3.0
  res3@xyLineColor       = "midnightblue"
  plot(1) = gsn_csm_xy (wks,time,TS_range,res2) ; draw plot  
  overlay(plot(1), gsn_csm_xy (wks,TS_time,TS_rbck,res))
  overlay(plot(1), gsn_csm_xy (wks,time,TS_diff,res3)) ; draw  plot
  
  res@trYMaxF = 0.05
  res@trYMinF = -0.32
  res@tmYLPrecision = 1
  res@tiYAxisString = "Precipitation (mm d~S~-1~N~)"

  res2=res
  res2@xyLineThicknessF = 0.5
  res2@xyMonoDashPattern = True
  res2@xyLineColor       = "lightblue"
  res2@gsnXYFillColors = (/"lightblue"/)
  res3=res
  res3@xyLineThicknessF = 3.0
  res3@xyLineColor       = "midnightblue"
  plot(2) = gsn_csm_xy (wks,time,PRECT_range,res2) ; draw plot
  overlay(plot(2), gsn_csm_xy (wks,PRECT_time,PRECT_rbck,res)) ; draw plot
  overlay(plot(2), gsn_csm_xy (wks,time,PRECT_diff,res3)) ; draw plot
  
  resP@tiYAxisFontHeightF=0.015
  newplot = gsn_attach_plots(baseplot,plot,res,resP)
  draw(baseplot)
  
  txres = True
  txres@txFontHeightF = 0.02
  txres@txFontColor="midnightblue"
  gsn_text_ndc(wks,"CESM1(WACCM)",0.73,0.88,txres)
  txres@txFontColor="red"
  gsn_text_ndc(wks,"GISS ModelE",0.73,0.84,txres)
    
  delete(wks)
  
  if (host.eq."modeling1") then
    muttcmd="mutt mikejmills@mac.com -a "+IMAGE_PATH+"/"+pdffile+".pdf -s "+pdffile+" < /dev/null"
    print(muttcmd)
    mutter=systemfunc(muttcmd)
    print(mutter)
  end if

  if (ostype.eq."darwin") then
    opencmd="open "+IMAGE_PATH+"/"+pdffile+".pdf"
    print(opencmd)
    open=systemfunc(opencmd)
    print(open)
  end if
  
;  print(PRECT_diff+" "+PRECT_range(0,:)+" "+PRECT_range(1,:))
