
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;========================================================================================================
; AVHRR and SAGE II data from Russell et al., 1996, Figure 1a, 19.5N:
; months since Jan 1, 1991
month_AVHRR195N = (/ 5.489,5.738,5.927,6.9,7.087,7.557,7.786,7.998,8.048,8.197,8.471,8.687, \
   8.946,9.143,9.41,9.592,9.844,10.08,10.309,10.541,10.59,10.766,11.039,11.45,11.663,11.855, \
   12.134,12.369,12.628,12.779,13.028,13.291,13.527,13.762,13.917,14.183,14.437,14.636,14.885, \
   15.082,15.36,15.551,15.795,16.122,16.225,16.434,16.518,16.693,16.961,17.187,17.403,17.864, \
   18.084,18.275,18.473,18.723,18.986,19.193,19.437,19.648,19.925,20.116,20.375,20.47,20.555, \
   20.736,21.004,21.149,21.511,21.755,21.955,22.157,22.406,22.64,22.799,23.12,23.305,23.603, \
   23.74,24.015,24.421,24.686,24.944,25.092,25.324,25.571,25.833,25.899,26.961,27.411,27.733, \
   29.416,29.697,30.143,30.371,30.421,30.605,31.102,31.276,31.532,31.798,32.03 /)
month_AVHRR195N@long_name = "Months after January 1, 1991"
AVHRR195N = (/  0,0.054,0.055,0.316,0.317,0.483,0.392,0.478,0.481,0.467,0.409,0.409,0.428,0.428, \
   0.28,0.279,0.26,0.314,0.298,0.332,0.329,0.297,0.315,0.26,0.3,0.315,0.281,0.294,0.277,0.224, \
   0.241,0.222,0.167,0.227,0.238,0.208,0.257,0.187,0.22,0.203,0.24,0.206,0.278,0.196,0.188,0.215, \
   0.216,0.168,0.168,0.186,0.186,0.222,0.277,0.242,0.242,0.224,0.258,0.205,0.201,0.168,0.207,0.222, \
   0.257,0.222,0.174,0.148,0.148,0.125,0.078,0.099,0.146,0.096,0.074,0.074,0.122,0.125,0.055,0.039, \
   0.054,0.019,0.054,0.057,0.075,0.021,0,0.033,0,0,0,0.033,0,0,0.091,0.022,0,0.005,0.053,0.049,0,0.055,0,0 /)
AVHRR195N@long_name="Stratospheric AOD in visible"

month_SAGE195Na = (/ 0.052,1.136,1.394,1.763,1.929,2.077,5.49,5.698,5.897,6.656,6.856 /)
SAGE195Na = (/ 0.00E+00,0.00E+00,0.018,0.019,0.014,0.00E+00,0.00E+00,0.032,0.037,0.038,0.149 /)

month_SAGE195Nb = (/9.207,9.488,9.669,9.895,10.785,11.097,11.443,11.593,12.187,12.305,12.473, \
   12.67,12.75,13.292,13.55,13.886,14.128,14.327,14.508,14.608,15.451,15.743,15.975,16.848, \
   16.98,17.141,17.441,17.622,19.21,19.39,20.511,20.808,22.115,22.434,22.884,23.062,23.575, \
   23.755,24.705,24.918,25.8,26.082,27.615,27.891,31.971 /)
SAGE195Nb = (/0.335,0.334,0.33,0.297,0.297,0.279,0.279,0.262,0.261,0.275,0.279,0.276,0.261,0.259, \
   0.241,0.242,0.274,0.28,0.274,0.262,0.259,0.278,0.279,0.279,0.274,0.245,0.238,0.203,0.204, \
   0.188,0.186,0.166,0.167,0.149,0.149,0.132,0.128,0.112,0.11,0.093,0.093,0.074,0.074,0.055,0.055 /)
;========================================================================================================

; AVHRR and SAGE II data from Russell et al., 1996, Figure 1b, 5S
; months since Jan 1, 1991
month_AVHRR5S = (/ 5.476,6.006,6.254,6.324,6.5,6.671,6.896,7.073,7.358,7.585,7.885,8.015,8.249,8.492, \
   8.682,8.994,9.26,9.371,9.724,10.113,10.319,10.557,11.029,11.276,11.651,11.9,12.115,12.322,12.631, \
   13.228,13.559,14.168,14.656,14.921,15.08,15.354,15.747,15.889,16.251,16.465,16.677,16.958,17.24, \
   17.586,18.041,18.31,18.56,18.764,19.062,19.23,19.466,19.689,19.904,20.127,20.569,20.643,20.841,21.05, \
   21.259,21.497,21.8,21.983,22.446,22.848,23.149,23.785,24.001,24.09,24.443,24.654,24.899,25.136,25.39, \
   25.558,25.853,26.978,27.253,27.373 /)
month_AVHRR5S@long_name = "Months after January 1, 1991"
AVHRR5S = (/ 0.00E+00,0.04,0.11,0.16,0.235,0.297,0.345,0.362,0.366,0.4,0.414,0.426,0.419,0.39,0.364,0.371, \
   0.355,0.338,0.319,0.297,0.279,0.244,0.233,0.232,0.205,0.207,0.214,0.215,0.196,0.159,0.159,0.185,0.15, \
   0.142,0.149,0.14,0.14,0.136,0.104,0.111,0.095,0.093,0.101,0.074,0.074,0.058,0.093,0.066,0.066,0.074, \
   0.066,0.065,0.083,0.083,0.048,0.048,0.056,0.037,0.037,0.046,0.028,0.035,0.019,0.018,0.028,0.028,0.042, \
   0.039,0.00E+00,0.00E+00,0.016,0.018,0.01,0.017,0.00E+00,0.00E+00,0.009,0.00E+00 /)
AVHRR5S@long_name="Stratospheric AOD in visible"

month_SAGE5S = (/ 0.017,5.544,5.737,6.239,6.746,6.873,7.369,7.549,7.789,8.162,8.44,8.766,8.991,9.199,9.534,9.698,9.89, \
   10.216,10.426,10.541,10.79,10.958,11.165,11.463,11.679,12.153,12.396,12.688,12.822,13.343,13.524,14.48, \
   14.686,15.568,15.773,16.276,16.49,16.986,17.099,17.364,17.587,17.844,18.713,19.038,19.245,19.44,19.483, \
   19.894,20.088,21.007,21.246,22.818,23.172,25.424,25.607,26.992,27.189,30.213,30.378,31.993 /)
SAGE5S = (/ 0.008,0.009,0.128,0.14,0.143,0.149,0.151,0.186,0.187,0.195,0.194,0.214,0.215,0.244,0.243,0.238, \
   0.205,0.205,0.208,0.214,0.214,0.222,0.223,0.224,0.179,0.178,0.188,0.188,0.195,0.196,0.178,0.176,0.167, \
   0.167,0.158,0.159,0.167,0.169,0.161,0.159,0.133,0.131,0.13,0.121,0.121,0.098,0.092,0.093,0.083,0.083, \
   0.074,0.074,0.066,0.065,0.056,0.057,0.047,0.047,0.038,0.038 /)
;========================================================================================================

; Data from Ansmann et al., 1997, Figure 4, 50-55N:
; days since eruption (Jun 15, 1991)
day_AVHRR50_55N = (/ 0.0,31.16,61.114,91.799,123.486,249.867,278.666,311.422,342.092,371.839,435.077,462.829, \
   497.326,649.978,712.347,736.33,805.42 /)
AVHRR50_55N = (/ 0.003,0.02,0.06,0.048,0.046,0.184,0.181,0.209,0.193,0.188,0.149,0.122,0.1,0.083,0.063,0.053,0.03 /)
day_AVHRR50_55N@long_name = "Days after eruption"
AVHRR50_55N@long_name = "Stratospheric AOD in visible"

day_SAGE_LIDAR50_55N = (/ 61.973,181.348,218.086,242.231,249.962,278.86,303.807,378.614,483.106,579.973,653.858, \
   737.363,802.638,928.393 /)
SAGE_LIDAR50_55N = (/ 0.037,0.101,0.161,0.195,0.207,0.23,0.226,0.202,0.153,0.115,0.094,0.077,0.064,0.038 /)

day_SAGE50_55N = (/ 0,32.094,62.886,89.817,121.581,146.638,155.314,159.203,164.033,168.896,189.267,207.579, \
   221.059,271.913,289.163,308.323,331.339,365.863,404.222,434.943,469.519,494.517,527.155,556.889,587.675, \
   624.112,635.602,649.013,682.608,716.162,742.069,776.64,802.548 /)
SAGE50_55N = (/ 0.01,0.013,0.026,0.034,0.049,0.069,0.077,0.088,0.095,0.109,0.159,0.175,0.184,0.173,0.164,0.152, \
   0.144,0.131,0.116,0.112,0.113,0.12,0.116,0.105,0.12,0.103,0.092,0.081,0.074,0.054,0.045,0.042,0.033 /)

day_LIDAR50_55N = (/ 62.023,92.737,123.544,154.5,186.188,217.04,248.896,279.738,310.501,341.043,373.555,403.302, \
   434.005,465.737,496.529,526.345,557.008,589.646,619.282,650.102,680.728,711.513,743.112,773.827,805.488,836.236, \
   866.975,897.73,928.436 /)
LIDARavg50_55N = (/ 0.048,0.043,0.059,0.112,0.11,0.14,0.181,0.21,0.218,0.17,0.133,0.127,0.119,0.129,0.146,0.16,0.14, \
   0.136,0.094,0.121,0.087,0.103,0.072,0.065,0.053,0.057,0.058,0.065,0.054 /)
LIDARmax50_55N = (/ 0.048,0.043,0.145,0.134,0.11,0.181,0.206,0.234,0.283,0.206,0.169,0.134,0.134,0.144,0.196,0.197, \
   0.186,0.204,0.121,0.175,0.087,0.144,0.09,0.102,0.065,0.076,0.058,0.065,0.054 /)
LIDARmin50_55N = (/ 0.048,0.043,0.031,0.088,0.11,0.092,0.141,0.154,0.18,0.122,0.095,0.094,0.103,0.09,0.113,0.123, \
   0.095,0.103,0.061,0.094,0.087,0.071,0.066,0.054,0.05,0.04,0.058,0.065,0.054 /)
;========================================================================================================
; Data from Aquila et al., JGR, 2012
month_global = (/ 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33 /)
AVHRR_global = (/ 0.013,0.017,0.097,0.167,0.181,0.16,0.164,0.16,0.144,0.128,0.138,0.132,0.127,0.126,\
                 0.122,0.1,0.097,0.067,0.052,0.054,0.044,0.032,0.023,0.026,0.027,0.024,0.028,0.013,0 /)
month_global@long_name="Months after January 1, 1991"
AVHRR_global@long_name="Stratospheric AOD in visible"

SAGE_global = (/ 0.007,0.012,0.046,0.092,0.107,0.12,0.133,0.134,0.139,0.145,0.134,0.126,0.125,0.112, \
                 0.11,0.108,0.092,0.088,0.079,0.069,0.069,0.061,0.054,0.049,0.045,0.038,0.039,0.037,0.032 /)
GEOS5_avg = (/ 0,0.04,0.151,0.19,0.202,0.199,0.188,0.171,0.156,0.143,0.136,0.128,0.118,0.109,0.099,0.089, \
               0.081,0.074,0.07,0.064,0.058,0.052,0.045,0.038,0.034,0.029,0.024,0.021,0.015 /)
GEOS5_max = (/ 0,0.04,0.162,0.2,0.205,0.203,0.193,0.184,0.174,0.168,0.164,0.155,0.146,0.136,0.126,0.116, \
               0.107,0.1,0.095,0.088,0.082,0.073,0.064,0.057,0.048,0.041,0.034,0.028,0.024 /)
GEOS5_min = (/ 0,0.04,0.141,0.182,0.196,0.196,0.178,0.153,0.13,0.114,0.104,0.096,0.082,0.07,0.058,0.052, \
               0.046,0.041,0.039,0.035,0.031,0.027,0.022,0.019,0.017,0.015,0.012,0.011,0.008 /)
 
;caseid = "f.e12.FSDW5.f19_f19.Pin10Tg.18-22km.001"
;caseid = "f.e12.FSDW5.f19_f19.Pin15Tg.20-27km.0-18N.96-118E.QBOon.001"
;caseid = "b.e12.B55TRW5CN.f19_g16.Pin10Tg.18-22km.QBOon.003"
;caseid = "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.003"
caseid = (/ \
 "b.e12.B55TRW5CN.f19_g16.Pin20Tg.16-18km.1col.QBOon.001", \
 "b.e12.B55TRW5CN.f19_g16.Pin20Tg.18-22km.1col.QBOon.001", \
 "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-24km.1col.QBOon.001", \
 "b.e12.B55TRW5CN.f19_g16.Pin20Tg.22-26km.1col.QBOon.001", \
 "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.1col.QBOon.001", \
 "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.001" \
; "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.002", \
; "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.003" \
; "b.e12.B55TRW5CN.f19_g16.Pin15Tg.20-27km.0-18N.96-118E.QBOon.avg", \
; "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-27km.0-18N.96-118E.QBOon.avg" \
; "f.e12.FSDW5.f19_f19.Pin15Tg.18-22km.001", \
; "f.e12.FSDW5.f19_f19.Pin15Tg.20-27km.0-18N.96-118E.QBOon.001" \
; "b.e12.B55TRW5CN.f19_g16.Pin10Tg.16-20km.QBOon.avg", \
; "b.e12.B55TRW5CN.f19_g16.Pin20Tg.16-18km.1col.QBOon.avg" \
  /)
  
rootpath =    "/glade/scratch/mmills/archive/"

nRuns = dimsizes(caseid)

do i = 0,nRuns-1
  i1=i+1
  print(i1+" of "+nRuns+": "+caseid(i))
  h0zmpath=rootpath+caseid(i)+"/atm/proc/h0zm/"
  h0gmpath=rootpath+caseid(i)+"/atm/proc/h0gm/"

  gmfile=addfile(h0gmpath+caseid(i)+".cam.h0.globavg.nc","r")
  AODvstrt_gm=gmfile->AODvstrt
  print("AODvstrt_gm:"+AODvstrt_gm)

  print("historyNames = systemfunc(ls "+ h0zmpath + "*))")
  historyNames = systemfunc("ls " + h0zmpath + "*")
  print(historyNames)
  print("history = addfiles(historyNames, r)")
  history = addfiles(historyNames, "r")
  print("ListSetType(history, cat)")
  ListSetType(history, "cat")
  ;  print(historyNames)
  print("historyNames:"+dimsizes(historyNames))

  print ("reading time")
  time    = addfiles_GetVar(history, historyNames, "time")
  print ("reading date")
  date    = addfiles_GetVar(history, historyNames, "date")
  print(time+" "+date)
  print ("reading AODvstrt")

  AODvstrt    = addfiles_GetVar(history, historyNames, "AODvstrt")
  printVarSummary(AODvstrt)

  day0=ind(date.eq.19910701)
  day1=ind(date.eq.19930101)
  ;print(day0+" "+time(day0))
  
  if (i.eq.0) then
    lat=AODvstrt&lat
    nLat=dimsizes(lat)
    nmonths = day1-day0+1
    day = time(day0:day1) - time(day0)
    day@long_name = "Days after eruption"
    month = ispan(6-day0, nmonths-day0+5, 1)
    month@long_name = "Months after January 1, 1991"
    aod=new( (/ nRuns, nmonths, nLat /), "double")
    aod_gm=new( (/ nRuns, nmonths /), "double")
    ;print(lat)
  end if  

  aod(i,:,:)  = AODvstrt(day0:day1,:)
  aod_gm(i,:) = AODvstrt_gm(day0:day1)
  
  delete(gmfile)
  delete(historyNames)
  delete(history)
  delete(time)
  delete(date)
  delete(AODvstrt)
  delete(AODvstrt_gm)
  delete(h0zmpath)
  delete(h0gmpath)
  delete(day0)
  delete(day1)

end do

aod@long_name="Stratospheric AOD in visible"

;print(month)
;print(time+" "+date+" "+month+" "+day+" "+aod(:,{10})+" "+aod(:,{15})+" "+aod(:,{19.5}))
;print(month+" "+aod(:,{10})+" "+aod(:,{15})+" "+aod(:,{19.5}))

  ; Plot the data.
  pdfpath="$IMG_PATH/"
  a = systemfunc("mkdir -p "+pdfpath)
  pdffile=pdfpath+"ExtVsObs"
  print("Creating file "+pdffile+".pdf")
  wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file
;  gsn_define_colormap(wks, "rainbow+white+gray")  ; choose a colormap

  res                   = True
  res@gsnPaperOrientation = "landscape"
  res@gsnMaximize         = True        ; must include w/ Paper Orientation
  res@vpHeightF             = 0.7              ; change aspect ratio
  res@vpWidthF              = 1.0
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@xyLineThicknessF  = 4.0
  
  res@xyLineColor       = "red"
;  res@tiXaxisString     = "Months since Jan 1, 1991"
;  res@tiYaxisString     = "Stratospheric AOD, visible"
  res@tiMainString     = "19.5N"
  plot = gsn_csm_xy(wks, month_AVHRR195N, AVHRR195N, res)

  res@xyLineColor       = "blue"
  overlay(plot, gsn_csm_xy(wks, month_SAGE195Na, SAGE195Na, res))
  overlay(plot, gsn_csm_xy(wks, month_SAGE195Nb, SAGE195Nb, res))

;  res@xyLineColor       = "green"
;  overlay(plot, gsn_csm_xy(wks, month, aod(:,{14.2}), res))
  res@xyMonoLineColor   = False
  res@xyLineColors       = (/"black","black","black","black","black","green4","brown"/)
  res@xyDashPatterns = (/4,3,2,1,0,0/)
;  overlay(plot, gsn_csm_xy(wks, month, aod(:,{10.4}), res))
;  res@xyLineColor       = "darkgreen"
  overlay(plot, gsn_csm_xy(wks, month, aod(:,:,{19.89}), res))

  draw(plot)
  frame(wks)
  
  res@xyMonoLineColor   = True
  res@xyLineColor       = "red"
  res@tiMainString     = "5.0S"
  plot = gsn_csm_xy(wks, month_AVHRR5S, AVHRR5S, res)

  res@xyLineColor       = "blue"
  overlay(plot, gsn_csm_xy(wks, month_SAGE5S, SAGE5S, res))
  
  res@xyMonoLineColor   = False
;  res@xyLineColors       = (/"black","black","black","black","green4","brown"/)
;  res@xyLineColor       = "black"
  overlay(plot, gsn_csm_xy(wks, month, aod(:,:,{-4.7}), res))
;  res@xyLineColor       = "black"
;  overlay(plot, gsn_csm_xy(wks, month, aod(:,{-10.4}), res))
;  res@xyLineColor       = "orange"
;  overlay(plot, gsn_csm_xy(wks, month, aod(:,{-0.95}), res))

  draw(plot)
  frame(wks)
  
;  res@tiXaxisString     = "Days since eruption"
  res@tiMainString     = "50N"
  res@xyLineColor       = "black"
  plot = gsn_csm_xy(wks, day, aod(:,:,{50.2}), res)

  res@xyMonoLineColor   = True
  res@xyLineColor       = "green"
  overlay(plot, gsn_csm_xy(wks, day_SAGE_LIDAR50_55N, SAGE_LIDAR50_55N, res))
  
  res@xyLineColor       = "red"
  overlay(plot, gsn_csm_xy(wks, day_AVHRR50_55N, AVHRR50_55N, res))

  res@xyLineColor       = "blue"
  overlay(plot, gsn_csm_xy(wks, day_SAGE50_55N, SAGE50_55N, res))

  res@xyLineColor       = "purple"
  overlay(plot, gsn_csm_xy(wks, day_LIDAR50_55N, LIDARavg50_55N, res))
  
;  res@xyLineColor       = "black"
;  overlay(plot, gsn_csm_xy(wks, day, aod(:,:,{50.2}), res))
;  overlay(plot, gsn_csm_xy(wks, day, aod(:,{52.1}), res))
;  overlay(plot, gsn_csm_xy(wks, day, aod(:,{54}), res))
;  res@xyLineColor       = "brown"
;  overlay(plot, gsn_csm_xy(wks, day, aod(:,{55.9}), res))

  draw(plot)
  frame(wks)
  
  
  res@xyLineColor       = "red"
;  res@tiXaxisString     = "Months since Jan 1, 1991"
  res@tiMainString     = "Global Average"
  plot = gsn_csm_xy(wks, month_global, AVHRR_global, res)

  res@xyLineColor       = "blue"
  overlay(plot, gsn_csm_xy(wks, month_global, SAGE_global, res))
  
;  res@xyLineColor       = "green4"
;  overlay(plot, gsn_csm_xy(wks, month_global, GEOS5_avg, res))
;  res@xyLineThicknessF  = 1.0
;  overlay(plot, gsn_csm_xy(wks, month_global, GEOS5_max, res))
;  overlay(plot, gsn_csm_xy(wks, month_global, GEOS5_min, res))
  
;  res@xyLineColor       = "black"
  res@xyMonoLineColor   = False
;  res@xyLineColors       = (/"black","black","black","black","green4","brown"/)
;  res@xyLineThicknessF  = 3.0
  overlay(plot, gsn_csm_xy(wks, month, aod_gm, res))

  draw(plot)
  frame(wks)
  
  delete(wks)
