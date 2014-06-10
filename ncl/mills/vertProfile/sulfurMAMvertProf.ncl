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

  ; Open the histories
  history = addfile(rootpath+caseid+"/atm/proc/h0gm/"+caseid+".cam.h0.globavg.nc","r")
    
  z     = history->Z3
  z     = (/ z / 1000. /)
  
  date=history->date
  startTime=ind(date.eq.year*10000+201)
  stopTime=ind(date.eq.(year+1)*10000+101)

; GASES
  ; Read in the gases (mol/mol)
  h2so4 = history->H2SO4
  hso3  = history->HSO3
  ocs   = history->OCS
  s     = history->S
  so    = history->SO
  so2   = history->SO2
  so3   = history->SO3
  hso3  = history->HSO3
  dms   = history->DMS
  
  ; Read in the sulfate (kg/kg)
  so4   = history->so4
  so4   = (/ so4*MWDAIR/MWNH4HSO4 /) ; kg/kg -> mol/mol
  
  sox   = so2 ; get metadata
  sox   = (/ s + so + so2 + so3 + hso3 /)
  
  totalgas = h2so4 ; get metadata
  totalgas = (/ totalgas + hso3 + ocs + s + so + so2 + so3 + dms /)
  
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
  R146      = history->R_146 ; rate constant 146: SO2 + OH + M ->  HSO3 + M 
  arhomass  = history->arhomass ; air density (g/cm3)
  arho = arhomass*AVOGAD/MWDAIR  ; air density (molec/cm3)
  arho@units="molec/cm3"
  
  r_jso3    = r_jso3/arho       ; reaction rate (molec/molec/s)
  r_jso3@units="molec/molec/s"
  r_SO2_OH    = r_SO2_OH/arho       ; reaction rate (molec/molec/s)
  r_SO2_OH@units="molec/molec/s"
  r_SO3_H2O    = r_SO3_H2O/arho       ; reaction rate (molec/molec/s)
  r_SO3_H2O@units="molec/molec/s"
  
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
  gso2 = dim_avg_Wrap(so2(lev|:, time|startTime:stopTime))
  gso2 = (/ gso2*1e12 /)
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
  R146_tavg =dim_avg_Wrap(R146(lev|:, time|startTime:stopTime))
  
  tau_SOx_tavg =dim_avg_Wrap(tau_SOx(lev|:, time|startTime:stopTime))
  tau_SOx_wr_tavg =dim_avg_Wrap(tau_SOx_wr(lev|:, time|startTime:stopTime))
  
  ; Open a workstation
  datestr = systemfunc("date +%y%m%d_%H%M")
  wks = gsn_open_wks("pdf", "$IMAGE_PATH/sulfurMAMvertProf_vmr_" + caseid +"_"+year+"avg" ); + "_c" + datestr)

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

  plot = new(4, graphic)
  
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
  dum=gsn_add_text(wks,plot(1),"SO2+OH+M->HSO3+M",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="lightblue"
  dum=gsn_add_text(wks,plot(1),"SO3+H2O->H2SO4",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="red"
  dum=gsn_add_text(wks,plot(1),"SO3+hv->SO2+O",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="blue"
  dum=gsn_add_text(wks,plot(1),"SO2 wet removal",xpos,ypos,txres)


  resxy@tiXAxisString         = "SOx Lifetime (days)"
  delete(resxy@trXMinF)            
  delete(resxy@trXMaxF)            
  resxy@trXMaxF               = 1.e6
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  plot(2) = gsn_csm_xy(wks, tau_SOx_tavg, gz, resxy)
  resxy@xyLineThicknesses     = (/ 4.0 /)
  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(2), gsn_csm_xy(wks, tau_SOx_wr_tavg, gz, resxy))
  
  xpos=2000.
  ypos=130.0
  ydif=8.0
  txres@txFontHeightF = 0.025                     ; font smaller. default big
  txres@txFontColor="black"
  dum=gsn_add_text(wks,plot(2),"SO3 oxidation",xpos,ypos,txres)
  ypos=ypos-ydif
  txres@txFontColor="blue"
  dum=gsn_add_text(wks,plot(2),"SO2 wet removal",xpos,ypos,txres)

  resxy@tiXAxisString         = "Rate constant ("+R146@units+")"
  delete(resxy@trXMinF)            
  delete(resxy@trXMaxF)            
;  resxy@trXMaxF               = 1.e6
  resxy@xyLineThicknesses     = (/ 8.0 /)
  resxy@xyLineColors          = (/ "black" /)
  plot(3) = gsn_csm_xy(wks, R146_tavg, gz, resxy)
  xpos=3.e-32
  ypos=60.0
  ydif=8.0
  txres@txFontColor="black"
  txres@txFontHeightF = 0.025                     ; font smaller. default big
  dum=gsn_add_text(wks,plot(3),"SO2+OH+M->HSO3+M",xpos,ypos,txres)

  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnPaperOrientation   = "portrait"
  gsn_panel(wks, plot, (/ 2, 2 /), resP)
;  draw(plot)
;  frame(wks)
  delete(wks)
  
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
  
caseid = getenv("caseid")

if (ismissing(caseid)) then
  print("set environment variable caseid before running")
else
  print("Using caseid: " + caseid)
  print("")
  
  year=stringtointeger(getenv("year"))
  if (ismissing(year)) then
    year=1990
  end if
  print("Averaging year: "+year)

  plotS(rootpath, caseid, year)
end if
