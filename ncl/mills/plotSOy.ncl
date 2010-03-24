load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "ncl-scripts/getHistoryPath.ncl"

undef("plotSOy")
procedure plotSOy(historyRun:string,pattern:string,day:integer)

local historyNames

begin

  historyPath = getHistoryPath(historyRun)
  
  historyPattern=historyPath+"allmonthly/*.cam2."+pattern+"*.nc"
;  sizePattern=historyRun+"/all.*.cam2.h0."+year+"*.nc"
  historyNames = systemfunc("ls " + historyPattern)
;  sizeNames = systemfunc("ls " + sizePattern)
  print(historyNames)
;  print(sizeNames)
  history = addfiles(historyNames, "r")
;  sizeHistory = addfiles(sizeNames, "r")
  ListSetType(history, "cat")
;  ListSetType(sizeHistory, "cat")
  SO2=addfiles_GetVar(history, historyNames, "SO2")
; SO2[time | 55] x [lev | 66] x [lat | 46] x [lon | 72]
  H2SO4=addfiles_GetVar(history, historyNames, "H2SO4")
  OCS=addfiles_GetVar(history, historyNames, "OCS")
  S=addfiles_GetVar(history, historyNames, "S")
  SO=addfiles_GetVar(history, historyNames, "SO")
  SO3=addfiles_GetVar(history, historyNames, "SO3")
  HSO3=addfiles_GetVar(history, historyNames, "HSO3")
  time = addfiles_GetVar(history, historyNames, "time")
  print(time)
  SOY=addfiles_GetVar(history, historyNames, "SOY")
;  Z3=addfiles_GetVar(history, historyNames, "Z3")
  
;  if (isfilevar(history,"PURSULMD")) then
;    PURSULVMR=addfiles_GetVar(history, historyNames, "PURSULVMR")
;  else
;    PURSULVMR=HSO3
;    PURSULVMR(:,:,:,:)=0.0
;  end if
  
  MIXSULVMR=addfiles_GetVar(history, historyNames, "MIXSULVMR")
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  nlev=dimsizes(lev)
  
;************************************************
; plotting parameters
;************************************************
  wks   = gsn_open_wks ("pdf","SOy_30N_"+historyRun+"_"+pattern)      ; open workstation

  res                   = True                     ; plot mods desired
  res@gsnPaperOrientation = "landscape"
  res@gsnMaximize       = True
  res@gsnBoxMargin      = 0.0
;  res@tiMainString      = "SO2 (pptv)"             ; add title
  res@trYReverse        = True                     ; reverse Y-axis
  res@trYLog            = True                     ; log Y-axis
  res@trXLog            = True                     ; log X-axis
;  res@trYMinF           = 0.3
;  res@trYMaxF           = 30.0
  res@trYMinF           = 5.0e-6
  res@trYMaxF           = 1013.0
  res@trXMinF           = 1.0
  res@trXMaxF           = 1000.0
  res@tiYAxisString     ="Pressure (hPa)"
  res@tiXAxisString     ="Volume mixing ratio (pptv)"
  
  res@tmYLMode   = "Explicit" ; Define own tick mark labels.
  res@tmYLValues = (/1000.,  100.,  10., 1., 0.1, 0.01, 0.001, 0.0001/)
  res@tmYLLabels = (/"1000","100", "10", "1", "0.1", "0.01", "0.001", "0.0001"/)
  res@tmYLMinorOn= False        ; No minor tick marks.
  
; If the maximium height is greater than 85 km, then theheight estimates from
; gsn_geop_hgt are not accurate. In this case, use a height axis which is just
; a table of values every 10km from 0 - 200km based on the US Standard
; Atmosphere.

  res@tmYRMode     = "Explicit"
  res@tmYUseLeft   = False
  res@tmYROn       = True
  res@tmYRLabelsOn = True
  res@tmYRValues   =  (/ 1013.25, 264.99, 55.293, 11.970, 2.8714, 0.79779, \
               0.21958, 5.2209e-2, 1.0524e-2, 1.8359e-3, 3.2011e-4, \
               7.1042e-5, 2.5382e-5, 1.2505e-5, 7.2028e-6, 4.5422e-6, \
               3.0395e-6, 2.1210e-6, 1.5271e-6, 1.1266e-6, 8.4736e-7 /)
  res@tmYRLabels   =  (/ "0", "10", "20", "30", "40", "50", \
               "60", "70", "80", "90", "100", \
               "110", "120", "130", "140", "150", \
               "160", "170", "180", "190", "200" /)

  res@xyLineThicknessF = 4.0

  data      = new((/6,nlev/),float)

  data(0,:) = dim_avg(SO2({day},:,{30},:))*1e12 ; longitudinal average
  print("SO2:"+min(data(0,:)))
  data(1,:) = dim_avg(H2SO4({day},:,{30},:))*1e12
  print("H2SO4:"+min(data(1,:)))
  data(2,:) = dim_avg(OCS({day},:,{30},:))*1e12
  print("OCS:"+min(data(2,:)))
  data(3,:) = dim_avg(SO3({day},:,{30},:))*1e12
  print("SO3:"+min(data(3,:)))
  data(4,:) = dim_avg(SOY({day},:,{30},:))*1e12
  print("SOY:"+min(data(4,:)))
;  data(5,:) = dim_avg(PURSULVMR({day},:,{30},:))*1e12
  data(5,:) = dim_avg(MIXSULVMR({day},:,{30},:))*1e12
  print("MIXSULVMR:"+min(data(5,:)))
  data = where (data.gt.0, data, data@_FillValue)
  print("DATA:"+min(data))

  res@xyDashPattern = 0               ; Make curves all solid
  res@xyLineColors = (/"Red","Blue","Green","Brown","Black","Orange" /)
  
;  z=dim_avg(Z3({day},:,{30},:))

  res@gsnDraw   = False
  res@gsnFrame  = False
  
  calldraw      = get_res_value(res,"gsnDraw", True)
  callframe     = get_res_value(res,"gsnFrame",True)
  maxbb         = get_bb_res(res)
  
  plot  = gsn_csm_xy (wks,data,lev,res) ; create plot
  
  rightaxis_string = create "right_axis" textItemClass wks
    "txString"      : "Height (km)"
    "txAngleF"      : 90.
  end create

  anno = NhlAddAnnotation(plot,rightaxis_string)
  setvalues anno
    "amZone"          : 3      ; Just outside plot area
    "amJust"          : "centercenter"
    "amSide"          : "right" 
    "amParallelPosF"  : 0.5
    "amOrthogonalPosF": 0.03
    "amResizeNotify"  : True     ; Resize if plot resized.
  end setvalues
  
  draw_and_frame(wks,plot,calldraw,callframe,0,maxbb)

end

;plotSOy("nongeo2","h0.1996-05",151+365)
;plotSOy("geo2","h0.1996-05",151+365)

;plotSOy("coag04","1995-01-01",1)
;plotSOy("coag04-2","1995-12-31",365+1)
;plotSOy("coag04-2","1996-12-26",728)
;plotSOy("nopur2_ncar","1995-07-05",185)

;plotSOy("125smin","monthly/*.cam2.h0.1997-11",244)
plotSOy("125smin","h0.1997-05",61)
