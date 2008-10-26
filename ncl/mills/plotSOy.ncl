load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotSOy")
procedure plotSOy(historyRun:string,year:string,day:integer)

local historyNames

begin

  historyPattern=historyRun+"/*.cam2."+year+"*.nc"
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
  wks   = gsn_open_wks ("pdf","SO2_30N_"+historyRun+"_d"+day)      ; open workstation

  res                   = True                     ; plot mods desired
  res@tiMainString      = "SO2 (pptv)"             ; add title
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
  res@xyLineThickness = 4.0

  data2      = new((/6,nlev/),float)

  data2(0,:) = dim_avg(SO2({day},:,{30},:))*1e12 ; longitudinal average
  print("SO2:"+min(data2(0,:)))
  data2(1,:) = dim_avg(H2SO4({day},:,{30},:))*1e12
  print("H2SO4:"+min(data2(1,:)))
  data2(2,:) = dim_avg(OCS({day},:,{30},:))*1e12
  print("OCS:"+min(data2(2,:)))
  data2(3,:) = dim_avg(SO3({day},:,{30},:))*1e12
  print("SO3:"+min(data2(3,:)))
  data2(4,:) = dim_avg(SOY({day},:,{30},:))*1e12
  print("SOY:"+min(data2(4,:)))
;  data2(5,:) = dim_avg(PURSULVMR({day},:,{30},:))*1e12
  data2(5,:) = dim_avg(MIXSULVMR({day},:,{30},:))*1e12
  print("MIXSULVMR:"+min(data2(5,:)))
  data = where (data.gt.0, data, data@_FillValue)
  print("DATA:"+min(data2))

  res@xyDashPattern = 0               ; Make curves all solid
  res@xyLineColors = (/"Red","Blue","Green","Brown","Black","Orange" /)
  
;  z=dim_avg(Z3({day},:,{30},:))
  plot  = gsn_csm_xy (wks,data2,lev,res) ; create plot

end

plotSOy("nongeo2","h0.1996-05",151+365)
;plotSOy("geo2","h0.1996-05",151+365q)

;plotSOy("coag04","1995-01-01",1)
;plotSOy("coag04-2","1995-12-31",365+1)
;plotSOy("coag04-2","1996-12-26",728)
;plotSOy("nopur2_ncar","1995-07-05",185)
