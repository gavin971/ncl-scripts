load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("plotS")
procedure plotS(rootpath, caseid, year)

begin

  AVOGAD  = 6.02214e23              ; Avogadro's number ~ molecules/mole
  MWDAIR  = 28.966                  ; molecular weight dry air ~ g/mole
  MWNH4HSO4 = 115.108999            ; molecular weight for MAM so4 species (ammonium sulfate)
  
  ; Hopfner observations
  Hopfner_DJF_SO2 = (/ 20.804,16.138,20.898,24.58,24.58,24.086,24.086,25.908,28.456,31.877, \
  36.355,39.073,39.548,38.839,37.918,39.311,41.561,45.995,50.905,57.364,65.084,70.742, \
  74.618,78.766,78.766,80.202 /)
  Hopfner_DJF_alt = (/ 15.121,15.974,18.021,18.97,20.435,22.03,23.01,24.257,25.163,26.032, \
  27.076,27.895,29.203,30.373,31.921,33.015,34.182,35.265,36.382,37.576,38.81,40.271,41.797, \
  43.208,44.248,44.906 /)

  Hopfner_MAM_SO2 = (/ 24.038,21.217,20.588,21.62,23.802,24.767,25.198,24.748,23.586,23.165, \
  24.179,26.579,29.772,33.524,37.327,41.062,43.871,46.066,47.473,47.438,47.438,47.438,48.335, \
  50.984,54.8,59.207,64.742,70.318,74.85,80.446,83.904,87.907,91.68,93.352,93.352,93.352,93.352 /)
  Hopfner_MAM_alt = (/ 15.062,15.721,16.517,17.343,18.116,18.945,19.848,20.674,21.55,22.376, \
  23.225,24.038,24.75,25.359,26.151,26.913,27.689,28.548,29.39,30.317,31.139,32.042,32.911, \
  33.727,34.553,35.339,36.142,36.974,37.763,38.536,39.428,40.211,41.147,42.02,42.892,43.745,44.621 /)

  Hopfner_JJA_SO2 = (/ 37.003,32.786,32.201,29.446,26.766,23.103,18.509,18.288,21.608,25.838, \
  32.614,39,39,39,37.62,38.074,40.891,46.007,50.652,59.532,65.785,71.94,77.668,81.985,85.509, \
  87.063 /)
  Hopfner_JJA_alt = (/ 15.145,16.024,17.332,18.499,19.646,20.488,21.95,22.879,23.966,24.788, \
  25.986,28.099,29.112,30.597,31.583,33.075,34.104,35.412,36.368,37.562,38.612,39.919,40.982, \
  42.333,43.477,44.911 /)

  Hopfner_SON_SO2 = (/ 27.806,24.127,22.855,22.855,23.536,23.963,23.236,21.011,19.666,17.863, \
  18.856,21.506,24.235,27.828,32.097,34.836,37.638,38.293,38.293,37.131,35.629,34.161,33.929, \
  34.574,37.778,40.816,46.343,51.562,57.499,62.544,67.225,71.765,75.754,78.538,79.486,80.931 /)
  Hopfner_SON_alt = (/ 15.14,15.849,16.695,17.507,18.199,18.647,19.49,20.389,21.138,21.887, \
  22.874,23.727,24.202,24.89,25.643,26.458,27.398,28.354,29.257,30.2,31.107,32.029,32.855, \
  33.915,34.781,35.47,36.453,37.138,37.938,38.841,39.64,40.6,41.533,42.402,43.355,44.311 /)

  Hopfner_DJF_std = (/ 7.355,9.309,8.786,9.105,10.485,9.303,10.11,10.658,9.414,9.857,11.704 \
  ,12.44,12.655,13.33,13.721,14.477,15.791,16.462,17.68,17.367,16.462,14.98,14.114,13.459, \
  14.38,15.249,16.66,17.68,18.639,19.315,20.226 /)
  Hopfner_DJF_std_alt = (/ 14.962,15.993,16.891,17.922,18.973,19.931,21.152,22.033,22.965, \
  24.042,24.933,26.018,27.082,27.957,29.337,30.212,31.156,32.363,32.975,34.037,35.17,36.255, \
  36.947,38.057,39.202,40.252,41.503,42.032,43.365,44.091,44.872 /)
  Hopfner_DJF_std2 = linint1 (Hopfner_DJF_std_alt, Hopfner_DJF_std, False, Hopfner_DJF_alt, 0)
  nalt=dimsizes(Hopfner_DJF_alt)
  Hopfner_DJF_std3 = conform_dims((/nalt,2/),Hopfner_DJF_std2,0)
  Hopfner_DJF_std3(:,0) = Hopfner_DJF_SO2 - Hopfner_DJF_std2
  Hopfner_DJF_std3(:,1) = Hopfner_DJF_SO2 + Hopfner_DJF_std2
  print(Hopfner_DJF_alt+" "+Hopfner_DJF_std3(:,0)+" "+Hopfner_DJF_std3(:,1))

  Hopfner_MAM_std = (/ 10.427,10.733,10.174,11.19,12.688,13.647,12.812,11.667,10.743,10.367, \
  11.502,11.299,10.912,10.774,11.025,11.885,13.455,14.765,16.117,16.815,16.815,16.617,16.323, \
  15.633,15.449,15.175,14.918,15.187,16.106,16.69,17.374,18.452,18.896,19.568,20.173,20.521, \
  21.234 /)
  Hopfner_MAM_std_alt = (/ 15.119,15.887,16.831,17.466,18.178,18.966,19.815,20.56,21.218, \
  22.119,22.911,23.719,24.564,25.435,26.393,27.204,27.79,28.574,29.423,30.188,31.069,31.944, \
  32.895,33.704,34.598,35.439,36.341,37.189,38.064,38.902,39.76,40.615,41.459,42.301,43.199, \
  44.02,44.952 /)
  Hopfner_MAM_std2 = linint1 (Hopfner_MAM_std_alt, Hopfner_MAM_std, False, Hopfner_MAM_alt, 0)
  nalt=dimsizes(Hopfner_MAM_alt)
  Hopfner_MAM_std3 = conform_dims((/nalt,2/),Hopfner_MAM_std2,0)
  Hopfner_MAM_std3(:,0) = Hopfner_MAM_SO2 - Hopfner_MAM_std2
  Hopfner_MAM_std3(:,1) = Hopfner_MAM_SO2 + Hopfner_MAM_std2
  print(Hopfner_MAM_alt+" "+Hopfner_MAM_std3(:,0)+" "+Hopfner_MAM_std3(:,1))

  Hopfner_JJA_std = (/ 9.222,11.262,13.372,14.013,13.283,12.006,9.955,10.134,11.08,11.707, \
  10.348,10.974,12.186,12.479,13.402,14.287,14.994,16.008,16.577,16.663,16.356,16.663,17.358, \
  18.203,20.594,22.838,25.456,27.319 /)
  Hopfner_JJA_std_alt = (/ 14.967,15.585,16.685,16.918,17.707,18.973,19.871,21.199,21.964, \
  23.926,25.822,27.119,28.382,28.938,29.463,30.817,31.984,33.282,34.373,35.726,36.834, \
  38.237,38.923,39.292,40.555,41.586,42.9,43.958 /)
  Hopfner_JJA_std2 = linint1 (Hopfner_JJA_std_alt, Hopfner_JJA_std, False, Hopfner_JJA_alt, 0)
  nalt=dimsizes(Hopfner_JJA_alt)
  Hopfner_JJA_std3 = conform_dims((/nalt,2/),Hopfner_JJA_std2,0)
  Hopfner_JJA_std3(:,0) = Hopfner_JJA_SO2 - Hopfner_JJA_std2
  Hopfner_JJA_std3(:,1) = Hopfner_JJA_SO2 + Hopfner_JJA_std2
  print(Hopfner_JJA_alt+" "+Hopfner_JJA_std3(:,0)+" "+Hopfner_JJA_std3(:,1))

  Hopfner_SON_std = (/ 6.921,8.605,9.619,8.211,7.334,8.018,8.396,7.334,6.363,7.772,9.104, \
  10.186,9.895,9.057,10.141,11.567,13.341,14.586,15.653,16.527,17.422,18.041,19.491,19.506, \
  18.614,18.069,17.126,15.572,14.086,13.441,13.204,13.441,14.002,14.511,15.759,17.422 /)
  Hopfner_SON_std_alt = (/ 15.019,15.624,16.02,16.845,17.906,19.093,20.021,21.128,22.026, \
  23.276,23.958,25.099,25.677,27.058,27.51,28.219,29.007,29.802,30.63,31.558,32.509,33.381, \
  34.262,35.117,36.224,37.073,37.944,38.716,39.541,40.016,40.572,41.489,42.281,43.103,44.2,44.908 /)
  Hopfner_SON_std2 = linint1 (Hopfner_SON_std_alt, Hopfner_SON_std, False, Hopfner_SON_alt, 0)
  nalt=dimsizes(Hopfner_SON_alt)
  Hopfner_SON_std3 = conform_dims((/nalt,2/),Hopfner_SON_std2,0)
  Hopfner_SON_std3(:,0) = Hopfner_SON_SO2 - Hopfner_SON_std2
  Hopfner_SON_std3(:,1) = Hopfner_SON_SO2 + Hopfner_SON_std2
  print(Hopfner_SON_alt+" "+Hopfner_SON_std3(:,0)+" "+Hopfner_SON_std3(:,1))
  
  ; Output from UCKA model
  UKCA_alt = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,\
  31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50/)
  
  UKCA_SO2_DJF = (/ 172.469,142.065,135.799,132.998,110.064,57.9482,38.4325,42.0034,46.4062, \
  49.3669,49.4707,46.4031,40.9751,33.8484,25.044,15.057,7.5267,2.54409,1.47825,1.32698,1.8107, \
  2.69877,4.02698,5.65229,7.48853,9.38924,10.7539,11.3277,11.028,9.83769,7.90527,5.97518,4.29498, \
  3.00708,2.24687,1.71124,1.31806,1.04642,0.880836,0.768996,0.728946,0.740852,0.774883,0.805098, \
  0.823877,0.802899,1.17002,2.30003,4.85315,12.4178 /)
  
  UKCA_SO2_MAM = (/ 118.917,105.88,95.3692,86.8596,72.4364,40.8215,26.5042,28.732,33.3789,37.9488, \
  40.8497,40.7436,37.0743,30.0337,20.2231,11.0847,5.15593,2.01584,1.50308,1.72093,2.40638,3.41914, \
  4.8586,6.61717,8.63871,10.8405,12.7628,14.0007,14.3977,13.8483,12.1959,9.77334,7.17865,4.70586, \
  2.82208,1.78274,1.16729,0.88523,0.793362,0.727227,0.681354,0.657198,0.667832,0.7506,1.09409, \
  1.64188,3.01318,4.76912,7.92577,11.8781 /)
  
  UKCA_SO2_JJA = (/ 110.697,88.7137,77.4081,61.1329,49.7367,34.2094,27.2087,32.0031,38.8663,45.7669, \
  51.3313,54.0317,52.3154,45.0092,31.1216,16.4288,6.6961,1.60263,1.07312,1.59298,2.58676,3.7749, \
  5.25426,7.02066,9.00204,11.1185,13.0769,14.3826,14.8612,14.3796,12.7186,10.139,7.42982,4.89351, \
  2.95416,1.82223,1.08589,0.676321,0.503633,0.375232,0.286291,0.224389,0.172615,0.149411,0.161901, \
  0.225629,0.344655,0.530392,0.845792,1.41572 /)
  
  UKCA_SO2_SON = (/ 133.143,107.068,82.7763,72.5234,62.6672,41.2643,33.3672,38.3535,44.3473,49.1865, \
  51.1622,48.6607,41.5143,30.903,19.2193,10.08,4.26289,1.4018,1.16616,1.56849,2.38261,3.48807, \
  5.08698,7.11499,9.53803,12.3112,15.098,17.1977,18.2343,17.897,15.6624,12.0374,8.26115,4.86385, \
  2.54053,1.42906,0.80771,0.531264,0.416072,0.325112,0.250639,0.206911,0.179045,0.183424,0.217064, \
  0.274655,0.331735,0.393479,0.4895,0.63328 /)
  
  ; Open the histories
;  history = addfile(rootpath+caseid+"/atm/proc/h0gm/"+caseid+".cam.h0.globavg.nc","r")
  history = addfile(rootpath+caseid(0)+"/atm/proc/h0gm/"+caseid(0)+".cam.h0.0-20N.nc","r")
  history1 = addfile(rootpath+caseid(1)+"/atm/proc/h0gm/"+caseid(1)+".cam.h0.0-20N.nc","r")
  history2 = addfile(rootpath+caseid(2)+"/atm/proc/h0gm/"+caseid(2)+".cam.h0.0-20N.nc","r")
    
  z     = history->Z3
  z     = (/ z / 1000. /)
  
  date=history->date
  startTime=ind(date.eq.year(0)*10000+201)
  stopTime=ind(date.eq.(year(0)+1)*10000+101)
  
  date1=history1->date
  startTime1=ind(date1.eq.year(1)*10000+201)
  date2=history2->date
  startTime2=ind(date2.eq.year(2)*10000+201)

; GASES
  ; Read in the gases (mol/mol)
  h2so4 = history->H2SO4
  ocs   = history->OCS
  s     = history->S
  so    = history->SO
  so2   = history->SO2
  so2_1 = history1->SO2
  so2_2 = history2->SO2
  so3   = history->SO3
  dms   = history->DMS
  T     = history->T ; temperature
  

  ; Read in the sulfate (kg/kg)
  so4   = history->so4
  so4   = (/ so4*MWDAIR/MWNH4HSO4 /) ; kg/kg -> mol/mol
  
  sox   = so2 ; get metadata
  sox   = (/ s + so + so2 + so3 /)
  
  totalgas = h2so4 ; get metadata
  totalgas = (/ totalgas + ocs + s + so + so2 + so3 + dms /)
  
  ; Calculate the total in pptv
  totals = totalgas
  totals = (/ totalgas + so4 /)
  totals = (/ totals * 1e12 /)
  totals@units = "pptv"
  totals@long_name = "Total Sulfur"

; REACTION RATES
  r_jso3    = history->r_jso3   ; reaction rate (molec/cm3/s)
  r_SO2_OH  = history->r_SO2_OH
  r_SO3_H2O = history->r_SO3_H2O
  DTWR_SO2  = history->DTWR_SO2 ; wet removal Neu scheme tendency (mol/mol/s)
  r_SO2_OH      = history->r_SO2_OH ; reaction rate: SO2 + OH + M ->  SO3 + M (cm6/molec2/s)
  jh2so4    = history->jh2so4 ; photolysis rate (/s)
  
  arhomass  = history->arhomass ; air density (g/cm3)
  arho = arhomass*AVOGAD/MWDAIR  ; air density (molec/cm3)
  arho@units="molec/cm3"
  
  r_jso3          = r_jso3/arho       ; reaction rate (molec/molec/s)
  r_jso3@units    = "molec/molec/s"
  r_SO2_OH        = r_SO2_OH/arho       ; reaction rate (molec/molec/s)
  r_SO2_OH@units  = "molec/molec/s"
  r_SO3_H2O       = r_SO3_H2O/arho       ; reaction rate (molec/molec/s)
  r_SO3_H2O@units = "molec/molec/s"
  
  A1=T ; get metadata
  A1=(/ arho*3.30E-31*(300./T)^4.30 /)
  B1=T ; get metadata
  B1=(/ A1/1.6e-12 /)
  SO2_OH_v1=T ; get metadata
  SO2_OH_v1=(/ (A1/(1+B1))*0.6^(1./(1+log(B1)^2.)) /)
  print("SO2_OH_v1: min="+min(SO2_OH_v1)+" max="+max(SO2_OH_v1))
  
  A2=T ; get metadata
  A2=(/ arho*3.0E-31*(300./T)^3.30 /)
  B2=T ; get metadata
  B2=(/ A2/1.5e-12 /)
  SO2_OH_v2=T ; get metadata
  SO2_OH_v2=(/ (A2/(1+B2))*0.6^(1./(1+log(B2)^2.)) /)
  print("SO2_OH_v2: min="+min(SO2_OH_v2)+" max="+max(SO2_OH_v2))
  
  tau_SOx = sox ; get metadata
  tau_SOx@_FillValue = default_fillvalue(typeof(tau_SOx))
  tau_SOx_wr=tau_SOx
  
  r_SO3_H2O@_FillValue = default_fillvalue(typeof(r_SO3_H2O))
  DTWR_SO2@_FillValue = default_fillvalue(typeof(DTWR_SO2))
  r_SO3_H2O = where(r_SO3_H2O.ne.0,r_SO3_H2O,r_SO3_H2O@_FillValue)
  DTWR_SO2 = where(DTWR_SO2.lt.0,-DTWR_SO2,DTWR_SO2@_FillValue) ; flip the sign
  
  tau_SOx = tau_SOx/r_SO3_H2O  ; SOx lifetime (s)
  tau_SOx = tau_SOx/3600./24   ; SOx lifetime (s->days)
  tau_SOx@units="days"
  
  tau_SOx_wr = tau_SOx_wr/DTWR_SO2  ; SOx lifetime (s)
  tau_SOx_wr = tau_SOx_wr/3600./24   ; SOx lifetime (s->days)
  tau_SOx_wr@units="days"

  
  ; Make a zonal and time average.
  gts  = dim_avg_Wrap(totals(lev|:, time|startTime:stopTime))
  gz   =dim_avg_Wrap(z(lev|:, time|startTime:stopTime))
  
  gocs = dim_avg_Wrap(ocs(lev|:, time|startTime:stopTime))
  gocs = (/ gocs*1e12 /)
  
  so2 = (/ so2*1e12 /)
  gso2 = dim_avg_Wrap(so2(lev|:, time|startTime:stopTime))
  gso2_djf = dim_avg_Wrap(so2(lev|:, time|startTime-1:startTime+1))
  gso2_mam = dim_avg_Wrap(so2(lev|:, time|startTime+2:startTime+4))
  gso2_jja = dim_avg_Wrap(so2(lev|:, time|startTime+5:startTime+7))
  gso2_son = dim_avg_Wrap(so2(lev|:, time|startTime+8:startTime+10))

  so2_1 = (/ so2_1*1e12 /)
  gso2_djf1 = dim_avg_Wrap(so2_1(lev|:, time|startTime1-1:startTime1+1))
  gso2_mam1 = dim_avg_Wrap(so2_1(lev|:, time|startTime1+2:startTime1+4))
  gso2_jja1 = dim_avg_Wrap(so2_1(lev|:, time|startTime1+5:startTime1+7))
  gso2_son1 = dim_avg_Wrap(so2_1(lev|:, time|startTime1+8:startTime1+10))

  so2_2 = (/ so2_2*1e12 /)
  gso2_djf2 = dim_avg_Wrap(so2_2(lev|:, time|startTime2-1:startTime2+1))
  gso2_mam2 = dim_avg_Wrap(so2_2(lev|:, time|startTime2+2:startTime2+4))
  gso2_jja2 = dim_avg_Wrap(so2_2(lev|:, time|startTime2+5:startTime2+7))
  gso2_son2 = dim_avg_Wrap(so2_2(lev|:, time|startTime2+8:startTime2+10))
  
  gsox = dim_avg_Wrap(sox(lev|:, time|startTime:stopTime))
  gsox = (/ gsox*1e12 /)
  gh2so4 = dim_avg_Wrap(h2so4(lev|:, time|startTime:stopTime))
  gh2so4 = (/ gh2so4*1e12 /)
  gdms = dim_avg_Wrap(dms(lev|:, time|startTime:stopTime))
  gdms = (/ gdms*1e12 /)
  gso4 = dim_avg_Wrap(so4(lev|:, time|startTime:stopTime))
  gso4 = (/ gso4*1e12 /)
  
  r_jso3_tavg = dim_avg_Wrap(r_jso3(lev|:, time|startTime:stopTime))
  r_SO2_OH_tavg = dim_avg_Wrap(r_SO2_OH(lev|:, time|startTime:stopTime))
  r_SO3_H2O_tavg = dim_avg_Wrap(r_SO3_H2O(lev|:, time|startTime:stopTime))
  DTWR_SO2_tavg =dim_avg_Wrap(DTWR_SO2(lev|:, time|startTime:stopTime))
  r_SO2_OH_tavg =dim_avg_Wrap(r_SO2_OH(lev|:, time|startTime:stopTime))
  SO2_OH_v1_tavg =dim_avg_Wrap(SO2_OH_v1(lev|:, time|startTime:stopTime))
  SO2_OH_v2_tavg =dim_avg_Wrap(SO2_OH_v2(lev|:, time|startTime:stopTime))
  jh2so4_tavg = dim_avg_Wrap(jh2so4(lev|:, time|startTime:stopTime))
  T_tavg = dim_avg_Wrap(T(lev|:, time|startTime:stopTime))
  
  tau_SOx_tavg =dim_avg_Wrap(tau_SOx(lev|:, time|startTime:stopTime))
  tau_SOx_wr_tavg =dim_avg_Wrap(tau_SOx_wr(lev|:, time|startTime:stopTime))
  
  ; Open a workstation
  datestr = systemfunc("date +%y%m%d_%H%M")
;  wks = gsn_open_wks("pdf", "$IMAGE_PATH/sulfurMAMvertProf_vmr_" + caseid +"_"+year+"avg" ); + "_c" + datestr)
  pdffile="$IMAGE_PATH/sulfurMAMvertProf_vmr_0-20N_"+year(0)+"_"+ caseid(0) +"_c" + datestr
  wks = gsn_open_wks("pdf", pdffile)

  minP = 0.
  maxP = 1000.  

  resxy                       = True
  resxy@gsnDraw               = False
  resxy@gsnFrame              = False
  resxy@gsnMaximize           = True
  resxy@vpWidthF              = 0.6
  resxy@vpHeightF             = 0.8
  resxy@xyLineColors          = (/ "black" /)
  resxy@xyDashPatterns        = (/ 0 /)
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@tiXAxisString         = "Total Sulfur (pptv)"
  resxy@tiYAxisString         = "Altitude (km)"
  resxy@tiYAxisFontHeightF    = 0.02
  resxy@tiMainFontHeightF     = 0.03
  resxy@trXMinF               = 1.
  resxy@trXMaxF               = 2e3
  resxy@trYMaxF               = 140.
  resxy@trXLog                = True
  resxy@tmXBMinorPerMajor     = 8
  resxy@tmXTMinorPerMajor     = 8

  plot = new(6, graphic)
  
; PLOT GASES
  plot(0) = gsn_csm_xy(wks, gts, gz, resxy)

  resxy@xyLineThicknesses     = (/ 6.0 /)
  resxy@xyLineColors          = (/ "purple" /)
  overlay(plot(0), gsn_csm_xy(wks, gsox, gz, resxy))

  resxy@xyLineThicknesses     = (/ 4.0 /)
  resxy@xyLineColors          = (/ "green" /)
  overlay(plot(0), gsn_csm_xy(wks, gocs, gz, resxy))

  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(0), gsn_csm_xy(wks, gh2so4, gz, resxy))

  resxy@xyLineColors          = (/ "green" /)
  overlay(plot(0), gsn_csm_xy(wks, gocs, gz, resxy))

  resxy@xyLineColors          = (/ "brown" /)
  overlay(plot(0), gsn_csm_xy(wks, gdms, gz, resxy))

  resxy@xyLineColors          = (/ "orange" /)
  overlay(plot(0), gsn_csm_xy(wks, gso4, gz, resxy))

  resxy@xyLineColors          = (/ "red" /)
  overlay(plot(0), gsn_csm_xy(wks, gso2, gz, resxy))

  xpos=500.0
  ypos=130.0
  ydif=8.0
  txres = True
  txres@txFontColor="black"
  txres@txFontHeightF = 0.025                     ; font smaller. default big
  dum=gsn_add_text(wks,plot(0),"Total S",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="green"
  dum=gsn_add_text(wks,plot(0),"OCS",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="brown"
  dum=gsn_add_text(wks,plot(0),"DMS",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="red"
  dum=gsn_add_text(wks,plot(0),"SO2",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="purple"
  dum=gsn_add_text(wks,plot(0),"SOx",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="blue"
  dum=gsn_add_text(wks,plot(0),"H2SO4",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="orange"
  dum=gsn_add_text(wks,plot(0),"SO4",xpos,ypos,txres)

; PLOT REACTION RATES
  delete(resxy@trXMinF)            
  delete(resxy@trXMaxF)            
  resxy@trXMinF               = 1.e-27
  resxy@trXMaxF               = 2.e-15
  resxy@tiXAxisString         = "Reaction rate ("+r_SO2_OH@units+")"
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  plot(1) = gsn_csm_xy(wks, r_SO2_OH_tavg, gz, resxy)

  resxy@xyLineThicknesses     = (/ 4.0 /)
  resxy@xyLineColors          = (/ "lightblue" /)
  overlay(plot(1), gsn_csm_xy(wks, r_SO3_H2O_tavg, gz, resxy))
  resxy@xyLineColors          = (/ "red" /)
  overlay(plot(1), gsn_csm_xy(wks, r_jso3_tavg, gz, resxy))
  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(1), gsn_csm_xy(wks, DTWR_SO2_tavg, gz, resxy))

  xpos=1.e-21
  ypos=130.0
  ydif=8.0
  txres@txFontColor="black"
  txres@txFontHeightF = 0.025                     ; font smaller. default big
  dum=gsn_add_text(wks,plot(1),"SO2+OH+M->SO3+M",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="lightblue"
  dum=gsn_add_text(wks,plot(1),"SO3+H2O->H2SO4",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="red"
  dum=gsn_add_text(wks,plot(1),"SO3+hv->SO2+O",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="blue"
  dum=gsn_add_text(wks,plot(1),"SO2 wet removal",xpos,ypos,txres)

  resxy@tiXAxisString         = "T (K)"
  delete(resxy@trXMinF)            
  delete(resxy@trXMaxF)            
  resxy@trXMaxF               = 290.
  resxy@trXMinF               = 190.
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  resxy@trXLog                = False
  plot(2) = gsn_csm_xy(wks, T_tavg, gz, resxy)

  resxy@tiXAxisString         = "SOx Lifetime (days)"
  delete(resxy@trXMinF)            
  delete(resxy@trXMaxF)            
;  resxy@trYMaxF               = 85.
  resxy@trYMaxF               = 40.
;  resxy@trXMaxF               = 1.e6
  resxy@trXMaxF               = 1.e2
  resxy@trXLog                = True
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  plot(3) = gsn_csm_xy(wks, tau_SOx_tavg, gz, resxy)
  resxy@xyLineThicknesses     = (/ 4.0 /)
  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(3), gsn_csm_xy(wks, tau_SOx_wr_tavg, gz, resxy))
  
  xpos=60.
  ypos=65.0
  ydif=5.0
  txres@txFontHeightF = 0.025                     ; font smaller. default big
  txres@txFontColor="black"
;  dum=gsn_add_text(wks,plot(3),"SO3 oxidation",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="blue"
;  dum=gsn_add_text(wks,plot(3),"SO2 wet removal",xpos,ypos,txres)

  resxy@tiXAxisString         = "Rate constant ("+r_SO2_OH@units+")"
  delete(resxy@trXMinF)            
  delete(resxy@trXMaxF)            
  resxy@trYMaxF               = 85.
  resxy@trYMaxF               = 40.
  resxy@trXMaxF               = 2.e-12
  resxy@trXMinF               = 2.e-16
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  plot(4) = gsn_csm_xy(wks, r_SO2_OH_tavg, gz, resxy)
  resxy@xyLineThicknesses     = (/ 4.0 /)
  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(4), gsn_csm_xy(wks, SO2_OH_v1_tavg, gz, resxy))
  resxy@xyLineThicknesses     = (/ 4.0 /)
  resxy@xyLineColors          = (/ "red" /)
  overlay(plot(4), gsn_csm_xy(wks, SO2_OH_v2_tavg, gz, resxy))

  xpos=1.e-13
  ypos=75.0
  ydif=5.0
  txres@txFontColor="black"
  txres@txFontHeightF = 0.025                     ; font smaller. default big
  dum=gsn_add_text(wks,plot(4),"SO2+OH, model",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="blue"
  dum=gsn_add_text(wks,plot(4),"JPL11",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="red"
  dum=gsn_add_text(wks,plot(4),"MAM",xpos,ypos,txres)

  resxy@tiXAxisString         = "jH2SO4"
  delete(resxy@trXMinF)            
  delete(resxy@trXMaxF)            
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  plot(5) = gsn_csm_xy(wks, jh2so4_tavg, gz, resxy)

  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnPaperOrientation   = "landscape"
  gsn_panel(wks, plot, (/ 2, 3 /), resP)

  delete(plot)
  plot = new(4, graphic)
  
  delete(resxy)
  resxy                       = True
  resxy@gsnDraw               = False
  resxy@gsnFrame              = False
  resxy@gsnMaximize           = True
  resxy@vpWidthF              = 1.0
  resxy@vpHeightF             = 0.8
  resxy@xyDashPatterns        = (/ 0 /)
  resxy@tiXAxisString         = "SO2 (pptv)"
  resxy@tiYAxisString         = "Altitude (km)"
  resxy@tiYAxisFontHeightF    = 0.02
  resxy@tiMainFontHeightF     = 0.03
  resxy@trXMinF               = 0.6
  resxy@trXMaxF               = 120.
  resxy@trYMinF               = 15.
  resxy@trYMaxF               = 45.
  resxy@trXLog                = True
  resxy@tmXBMinorPerMajor     = 8
  resxy@tmXTMinorPerMajor     = 8
  
; PLOT SO2 DJF 
  resxy@tiMainString         = "DJF, 0-20N"
  resxy@xyDashPatterns        = (/ 1 /)
  resxy@xyLineThicknesses     = (/ 2.0 /)
  resxy@xyLineColors          = (/ "red" /)
  printVarSummary(Hopfner_DJF_std3)
  plot(0) = gsn_csm_xy(wks, Hopfner_DJF_std3(:,0), Hopfner_DJF_alt, resxy)
  overlay(plot(0), gsn_csm_xy(wks, Hopfner_DJF_std3(:,1), Hopfner_DJF_alt, resxy))

  delete(resxy@xyLineThicknesses)
  delete(resxy@xyLineColors)
  resxy@xyDashPatterns        = (/ 0 /)
  resxy@xyLineThicknesses     = (/ 3.0 /)
  resxy@xyLineColors          = (/ "purple" /)
;  overlay(plot(0), gsn_csm_xy(wks, UKCA_SO2_DJF, UKCA_alt, resxy))

  resxy@xyLineThicknesses     = (/ 5.0 /)
  resxy@xyLineColors          = (/ "red" /)
  overlay(plot(0), gsn_csm_xy(wks, Hopfner_DJF_SO2, Hopfner_DJF_alt, resxy))

  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  overlay(plot(0), gsn_csm_xy(wks, gso2_djf, gz, resxy))

  resxy@xyLineThicknesses     = (/ 3.0 /)
  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(0), gsn_csm_xy(wks, gso2_djf1, gz, resxy))
  resxy@xyLineColors          = (/ "green" /)
  overlay(plot(0), gsn_csm_xy(wks, gso2_djf2, gz, resxy))

  xpos=10.0
  ypos=42.0
  ydif=4.0
  txres = True
  txres@txFontColor="black"
  txres@txFontHeightF = 0.025                     ; font smaller. default big
  dum=gsn_add_text(wks,plot(0),"WACCM",xpos,ypos,txres)
  xpos=30.0
  ypos=40.0
  txres@txFontColor="blue"
;  dum=gsn_add_text(wks,plot(0),"QY=1",xpos,ypos,txres)
  txres@txFontColor="green"
  xpos=2.0
  ypos=41.0
;  dum=gsn_add_text(wks,plot(0),"old JH2SO4",xpos,ypos,txres)
  
  xpos=60.0
  ypos=20.0
  txres@txFontColor="red"
  dum=gsn_add_text(wks,plot(0),"MIPAS",xpos,ypos,txres)
  
  xpos=1.0
  ypos=35.0
  txres@txFontColor="purple"
;  dum=gsn_add_text(wks,plot(0),"UKCA",xpos,ypos,txres)

; PLOT SO2 MAM
  resxy@tiMainString         = "MAM, 0-20N"
  resxy@xyDashPatterns        = (/ 1 /)
  resxy@xyLineThicknesses     = (/ 2.0 /)
  resxy@xyLineColors          = (/ "red" /)
  printVarSummary(Hopfner_MAM_std3)
  plot(1) = gsn_csm_xy(wks, Hopfner_MAM_std3(:,0), Hopfner_MAM_alt, resxy)
  overlay(plot(1), gsn_csm_xy(wks, Hopfner_MAM_std3(:,1), Hopfner_MAM_alt, resxy))

  resxy@xyDashPatterns        = (/ 0 /)
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  overlay(plot(1), gsn_csm_xy(wks, gso2_mam, gz, resxy))

  resxy@xyLineThicknesses     = (/ 3.0 /)
  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(1), gsn_csm_xy(wks, gso2_mam1, gz, resxy))
  resxy@xyLineColors          = (/ "green" /)
  overlay(plot(1), gsn_csm_xy(wks, gso2_mam2, gz, resxy))

  resxy@xyLineThicknesses     = (/ 3.0 /)
  resxy@xyLineColors          = (/ "purple" /)
;  overlay(plot(1), gsn_csm_xy(wks, UKCA_SO2_MAM, UKCA_alt, resxy))

  resxy@xyLineThicknesses     = (/ 5.0 /)
  resxy@xyLineColors          = (/ "red" /)
  overlay(plot(1), gsn_csm_xy(wks, Hopfner_MAM_SO2, Hopfner_MAM_alt, resxy))

; PLOT SO2 JJA
  resxy@tiMainString         = "JJA, 0-20N"
  resxy@xyDashPatterns        = (/ 1 /)
  resxy@xyLineThicknesses     = (/ 2.0 /)
  resxy@xyLineColors          = (/ "red" /)
  printVarSummary(Hopfner_JJA_std3)
  plot(2) = gsn_csm_xy(wks, Hopfner_JJA_std3(:,0), Hopfner_JJA_alt, resxy)
  overlay(plot(2), gsn_csm_xy(wks, Hopfner_JJA_std3(:,1), Hopfner_JJA_alt, resxy))

  resxy@xyDashPatterns        = (/ 0 /)
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  overlay(plot(2), gsn_csm_xy(wks, gso2_jja, gz, resxy))

  resxy@xyLineThicknesses     = (/ 3.0 /)
  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(2), gsn_csm_xy(wks, gso2_jja1, gz, resxy))
  resxy@xyLineColors          = (/ "green" /)
  overlay(plot(2), gsn_csm_xy(wks, gso2_jja2, gz, resxy))

  resxy@xyLineThicknesses     = (/ 3.0 /)
  resxy@xyLineColors          = (/ "purple" /)
;  overlay(plot(2), gsn_csm_xy(wks, UKCA_SO2_JJA, UKCA_alt, resxy))

  resxy@xyLineThicknesses     = (/ 5.0 /)
  resxy@xyLineColors          = (/ "red" /)
  overlay(plot(2), gsn_csm_xy(wks, Hopfner_JJA_SO2, Hopfner_JJA_alt, resxy))

; PLOT SO2 SON
  resxy@tiMainString         = "SON, 0-20N"
  resxy@xyDashPatterns        = (/ 1 /)
  resxy@xyLineThicknesses     = (/ 2.0 /)
  resxy@xyLineColors          = (/ "red" /)
  printVarSummary(Hopfner_SON_std3)
  plot(3) = gsn_csm_xy(wks, Hopfner_SON_std3(:,0), Hopfner_SON_alt, resxy)
  overlay(plot(3), gsn_csm_xy(wks, Hopfner_SON_std3(:,1), Hopfner_SON_alt, resxy))

  resxy@xyDashPatterns        = (/ 0 /)
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  overlay(plot(3), gsn_csm_xy(wks, gso2_son, gz, resxy))

  resxy@xyLineThicknesses     = (/ 3.0 /)
  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(3), gsn_csm_xy(wks, gso2_son1, gz, resxy))
  resxy@xyLineColors          = (/ "green" /)
  overlay(plot(3), gsn_csm_xy(wks, gso2_son2, gz, resxy))

  resxy@xyLineThicknesses     = (/ 3.0 /)
  resxy@xyLineColors          = (/ "purple" /)
;  overlay(plot(3), gsn_csm_xy(wks, UKCA_SO2_SON, UKCA_alt, resxy))

  resxy@xyLineThicknesses     = (/ 5.0 /)
  resxy@xyLineColors          = (/ "red" /)
  overlay(plot(3), gsn_csm_xy(wks, Hopfner_SON_SO2, Hopfner_SON_alt, resxy))

  gsn_panel(wks, plot, (/ 2, 2 /), resP)
  
  delete(wks)
  print(pdffile+".pdf")
  
end

print("running sulfurMAMvertProf.ncl...")

host   = systemfunc("hostname")
print("host: "+host)

rootpath=getenv("rootpath")
if (ismissing(rootpath)) then
  if (host.eq."modeling1") then
    rootpath="/data2/mmills/"
  else
    rootpath = "/glade/scratch/mmills/archive/"
  end if
end if
print("Using rootpath: " + rootpath)
  
;caseid = getenv("caseid")
;caseid=(/"b.e12.B55TRW5CN.f19_g16.1989Tr.QBOon.NeuSO2.LyA.expJH2SO4.bkgd.2", \
;         "b.e12.B55TRW5CN.f19_g16.1989Tr.QBOon.NeuSO2.LyA.expJH2SO4.QY1.bkgd.1", \
;         "b.e12.B55TRW5CN.f19_g16.1989Tr.QBOon.NeuSO2.LyA.bkgd.1"/)
;year=(/1992,1992,1990/)
;caseid=(/"b.e13.B20TRW5TCN.f19_g16.beta13_stratmodal.004.expfix", \
;         "b.e12.B55TRW5CN.f19_g16.1989Tr.QBOon.NeuSO2.LyA.expJH2SO4.bkgd.2"/)
;year=(/1989,1992/)
caseid=(/"b.e13.B20TRW5TCN.f19_g16.beta13_stratmodal.004.expfix", \
         "b.e13.B20TRW5TCN.f19_g16.beta13_stratmodal.004.expfix", \
         "b.e13.B20TRW5TCN.f19_g16.beta13_stratmodal.004.expfix"/)
year=(/1988,1989,1990/)

;if (ismissing(caseid)) then
;  print("set environment variable caseid before running")
;else
  print("Using caseid: " + caseid + " year "+year)
  print("")
  
;  year=stringtointeger(getenv("year"))
;  if (ismissing(year)) then
;    year=1990
;  end if
;  print("Averaging year: "+year)

  plotS(rootpath, caseid, year)
;end if
