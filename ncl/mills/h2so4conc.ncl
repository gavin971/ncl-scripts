BOLTZ   = 1.38065e-16      ; Boltzmann's constant ~ erg/K/molecule

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "ncl-scripts/getHistoryPath.ncl"

undef("h2so4conc")
procedure h2so4conc(historyRun:string,prefix:string)

;  inputs: historyRun, i.e. "125smin"
;          prefix, i.e. "all.125smin.cam2.h3.1997-[01][09]"

local historyPath, historyPathD, historyNames, history

begin
  historyPath = getHistoryPath(historyRun)
  historyPathD = historyPath + "alldaily/"
  historyNames = systemfunc("ls " + historyPathD + prefix + "*")
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  print(historyNames)
  print("historyNames:"+dimsizes(historyNames))
  
  time    = addfiles_GetVar(history, historyNames, "time")
  nTime   = dimsizes(time)
  print("Number of Days:"+nTime)
  date = addfiles_GetVar(history, historyNames, "date")
  print("Date range:"+date(0)+" to "+date(nTime-1))

  lev     = history[0]->lev
  nLevs   = dimsizes(lev)
  lat     = history[0]->lat
  nLats   = dimsizes(lat)
  lon     = history[0]->lon
  nLons   = dimsizes(lon)
  
  ; Determine an approximate dz/dp.
  hyam 	= history[0]->hyam
  hybm 	= history[0]->hybm
  hyai 	= history[0]->hyai
  hybi 	= history[0]->hybi
  P0 	  = history[0]->P0
  PS 	  = history[0]->PS
  
  print("reading H2SO4")
  H2SO4    = addfiles_GetVar(history, historyNames, "H2SO4")
  print("reading H2SO4eq")
  H2SO4eq    = addfiles_GetVar(history, historyNames, "H2SO4EQ")
  print("reading P")
  P    = addfiles_GetVar(history, historyNames, "PRESSURE")
  print("reading T")
  T    = addfiles_GetVar(history, historyNames, "T")
  print("reading Z3")
  Z3    = addfiles_GetVar(history, historyNames, "Z3")

; Calculate number density of air molecules from ideal gas law
  print("calculating AIR")
  AIR  = P ; Get metadata
  AIR = AIR / BOLTZ / T
  print("calculating H2SO4 conc")
  H2SO4conc = H2SO4 ; Get metadata
  H2SO4eqconc = H2SO4eq ; Get metadata
  
  H2SO4conc = H2SO4conc * AIR 
  H2SO4eqconc = H2SO4eqconc * AIR 
  
  H2SO4ave = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(H2SO4conc(lev|:, time|:, {lat|38:46}, lon|:))))
  H2SO4eqave = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(H2SO4eqconc(lev|:, time|:, {lat|38:46}, lon|:))))
  Z        = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(Z3(lev|:, time|:, {lat|38:46}, lon|:))))
  Z = Z / 1000. ; convert m -> km
  
  H2SO4min = H2SO4ave ; Get metadata
  H2SO4max = H2SO4ave ; Get metadata
  H2SO4min = dim_min(dim_min(dim_min(H2SO4conc(lev|:, time|:, {lat|38:46}, lon|:))))
  H2SO4max = dim_max(dim_max(dim_max(H2SO4conc(lev|:, time|:, {lat|38:46}, lon|:))))
  
  H2SO4eqmin = H2SO4eqave ; Get metadata
  H2SO4eqmax = H2SO4eqave ; Get metadata
  H2SO4eqmin = dim_min(dim_min(dim_min(H2SO4eqconc(lev|:, time|:, {lat|38:46}, lon|:))))
  H2SO4eqmax = dim_max(dim_max(dim_max(H2SO4eqconc(lev|:, time|:, {lat|38:46}, lon|:))))
  
; when v 5.1.1 comes out...
;  H2SO4min = dim_min_n(H2SO4conc(time|:, {lat|38:46}, lon|:, lev|:),(/0,1,2/))
;  H2SO4max = dim_max_n(H2SO4conc(time|:, {lat|38:46}, lon|:, lev|:),(/0,1,2/))
;  H2SO4ave = dim_avg_n_Wrap(H2SO4conc(time|:, {lat|38:46}, lon|:, lev|:),(/0,1,2/))
;  Z        = dim_avg_n_Wrap(Z3(time|:, {lat|38:46}, lon|:, lev|:),(/0,1,2/))

  wks = gsn_open_wks("pdf", "H2SO4conc.h3."+historyRun)

  res = True
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
    
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "portrait"
;  res@gsnXYFillColors = "red" ; fill color to use between two adjacent curves in an XY plot (gsn_csm_xy), doesn't work for vertical plots
  
  res@xyMonoLineThickness = False
  res@xyLineThicknesses = (/2.0,4.0,2.0,2.0,4.0,2.0/)
  res@xyLineColors      = (/"red","brown","red","green","green4","green"/); change line color
  res@xyDashPatterns=(/0,0,0,0,0,0/)
  res@trXLog            = True                     ; log X-axis
  res@trYMaxF           = 100.0                    ; 100 km upper bound
;  res@xyComputeXMin     = True
  res@trXMinF           = 0.1
  res@trXMaxF           = 1e9
  res@tiYAxisString     ="Geopotential Height (km)"
  res@tiXAxisString     ="H~B~2~N~SO~B~4~N~ concentration (cm~S~-3~N~)"

;  print("lev:"+dimsizes(lev))
;  print("min:"+dimsizes(H2SO4min))
;  print("ave:"+dimsizes(H2SO4ave))
;  print("max:"+dimsizes(H2SO4max))
;  print("Z:"+dimsizes(Z))
  
  data= new((/6,dimsizes(lev)/),float)
  data(0,:)=H2SO4min(:)
  data(1,:)=H2SO4ave(:)
  data(2,:)=H2SO4max(:)
  data(3,:)=H2SO4eqmin(:)
  data(4,:)=H2SO4eqave(:)
  data(5,:)=H2SO4eqmax(:)
  
  plot = gsn_csm_xy(wks, data, Z, res)
  
;
; The fill_xy2 method was deprecated in V5.1.0. See instead the
; new resource "gsnXYFillColors".
;
;  plot  = fill_xy2(wks,plot(0),time,mnmx(2,:),mnmx(3,:),(/0.64,0.71,0.8/),\
;  (/0.64,0.71,0.8/))
;  plot  = fill_xy2(wks,plot(0),time,mnmx(0,:),mnmx(1,:),(/1.0,0.71,0.76/),\
;  (/1.0,0.71,0.76/))

  a = addfile(historyPath+"H2SO4conc.nc","c")
  a->lev=lev
  a->date=date
  a->Z=Z
  a->H2SO4min=H2SO4min
  a->H2SO4ave=H2SO4ave
  a->H2SO4max=H2SO4max
  a->H2SO4eqmin=H2SO4eqmin
  a->H2SO4eqave=H2SO4eqave
  a->H2SO4eqmax=H2SO4eqmax
  
end

;h2so4conc("125smin","all.125smin.cam2.h3.1997-[01][09]") 
;h2so4conc("nogeoPurSpin1","all.nogeoPurSpin1.cam2.h3.1995-[01][09]") 
h2so4conc("nogeo400SO2","all.nogeo400SO2.cam2.h3.1998-[01][09]") 
