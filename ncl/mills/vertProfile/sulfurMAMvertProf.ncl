load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("plotS")
procedure plotS(rootpath, caseid, year)

begin
  MWDAIR  = 28.966                  ; molecular weight dry air ~ kg/kmole
  MWNH4HSO4 = 115.108999            ; molecular weight for MAM so4 species (ammonium sulfate)

  ; Open the histories
  history = addfile(rootpath+caseid+"/atm/proc/h0gm/"+caseid+".cam.h0.globavg.nc","r")
    
  z     = history->Z3
  z     = (/ z / 1000. /)
  
  date=history->date
  startTime=ind(date.eq.year*10000+201)
  stopTime=ind(date.eq.(year+1)*10000+101)

  ; Read in the gases (mol/mol)
  h2so4 = history->H2SO4
  hso3  = history->HSO3
  ocs   = history->OCS
  s     = history->S
  so    = history->SO
  so2   = history->SO2
  so3   = history->SO3
  dms   = history->DMS
  
  ; Read in the sulfate (kg/kg)
  so4   = history->so4
  so4   = (/ so4*MWDAIR/MWNH4HSO4 /) ; kg/kg -> mol/mol
  
  totalgas = h2so4 ; get metadata
  totalgas = (/ totalgas + hso3 + ocs + s + so + so2 + so3 + dms /)
  
    
  ; Calculate the total in pptv
  totals = totalgas
  totals = (/ totalgas + so4 /)
  totals = (/ totals * 1e12 /)
  totals@units = "pptv"
  totals@long_name = "Total Sulfur"

  ; Make a zonal and time average.
  gts  = dim_avg_Wrap(totals(lev|:, time|startTime:stopTime))
  gz   =dim_avg_Wrap(z(lev|:, time|startTime:stopTime))
  
  gocs = dim_avg_Wrap(ocs(lev|:, time|startTime:stopTime))
  gocs = (/ gocs*1e12 /)
  gso2 = dim_avg_Wrap(so2(lev|:, time|startTime:stopTime))
  gso2 = (/ gso2*1e12 /)
  gso3 = dim_avg_Wrap(so3(lev|:, time|startTime:stopTime))
  gso3 = (/ gso3*1e12 /)
  gh2so4 = dim_avg_Wrap(h2so4(lev|:, time|startTime:stopTime))
  gh2so4 = (/ gh2so4*1e12 /)
  gdms = dim_avg_Wrap(dms(lev|:, time|startTime:stopTime))
  gdms = (/ gdms*1e12 /)
  gso4 = dim_avg_Wrap(so4(lev|:, time|startTime:stopTime))
  gso4 = (/ gso4*1e12 /)
   
  ; Open a workstation
  datestr = systemfunc("date +%y%m%d_%H%M")
  wks = gsn_open_wks("pdf", "$IMAGE_PATH/sulfurMAMvertProf_vmr_" + caseid +"_"+year+"avg" ); + "_c" + datestr)
;  gsn_define_colormap(wks,"BlWhRe")       
      
  resP                       = True
  resP@gsnMaximize           = True
;  resP@gsnPaperOrientation   = "portrait"
  resP@gsnPanelLabelBar     = False

  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@gsnMaximize           = True
;  res@vpWidthF              = 0.4
;  res@vpHeightF             = 0.8
      
;  res@tiMainFontHeightF     = 0.02
   

  res@gsnSpreadColors       = True 
  res@cnFillOn              = True
;  res@cnLevelSelectionMode  = "ManualLevels"
;  res@cnMinLevelValF        = -2.
;  res@cnMaxLevelValF        =  2.
;  res@cnLevelSpacingF       =  .2     
  res@cnLevelSelectionMode = "ExplicitLevels" 
  res@cnLevels             = (/ 1, 10, 20, 40, 80, 1e2, 2e2, 4e2, 8e2, 1e3, 2e3 /)
  res@lbLabelBarOn          = True
  res@cnInfoLabelOn         = False
  res@cnLineLabelsOn        = True
  res@cnLinesOn             = True

  res@lbTitleOn            = True                ; turn on title
  res@lbTitlePosition      = "Right"              ; title position
  res@lbTitleFontHeightF   = .025                ; make title smaller
  res@lbTitleDirection     = "Across"             ; title direction
  res@lbTitleString = "pptv"
  res@lbLabelStride        = 1
;  res@txString = "SPE - No SPE"

  minP = 0.
  maxP = 1000.  

  resxy                       = True
  resxy@gsnDraw               = False
  resxy@gsnFrame              = False
  resxy@gsnMaximize           = True
;  resxy@vpWidthF              = 0.8
;  resxy@vpHeightF             = 0.2
  resxy@xyLineColors          = (/ "black" /)
  resxy@xyDashPatterns        = (/ 0 /)
  resxy@xyLineThicknesses     = (/ 6.0 /)
  resxy@tiXAxisString         = "Total Sulfur (pptv)"
  resxy@tiYAxisString         = "Altitude (km)"
  resxy@tiYAxisFontHeightF    = 0.02
  resxy@tiMainFontHeightF     = 0.03
  resxy@trXMinF               = 1.
  resxy@trXMaxF               = 4e3
  resxy@trXLog                = True
  resxy@tmXBMinorPerMajor     = 8
  resxy@tmXTMinorPerMajor     = 8

  plot = gsn_csm_xy(wks, gts, gz, resxy)

  resxy@xyLineThicknesses     = (/ 4.0 /)
  resxy@xyLineColors          = (/ "green" /)
  overlay(plot, gsn_csm_xy(wks, gocs, gz, resxy))

  resxy@xyLineColors          = (/ "purple" /)
  overlay(plot, gsn_csm_xy(wks, gso3, gz, resxy))

  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot, gsn_csm_xy(wks, gh2so4, gz, resxy))

  resxy@xyLineColors          = (/ "green" /)
  overlay(plot, gsn_csm_xy(wks, gocs, gz, resxy))

  resxy@xyLineColors          = (/ "yellow" /)
  overlay(plot, gsn_csm_xy(wks, gdms, gz, resxy))

  resxy@xyLineColors          = (/ "orange" /)
  overlay(plot, gsn_csm_xy(wks, gso4, gz, resxy))

  resxy@xyLineColors          = (/ "red" /)
  overlay(plot, gsn_csm_xy(wks, gso2, gz, resxy))

  xpos=1000
  ypos=130
  ydif=10
  txres = True
  txres@txFontColor="black"
  txres@txFontHeightF = 0.03                     ; font smaller. default big
  gsn_text(wks,plot,"Total S",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="green"
  gsn_text(wks,plot,"OCS",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="yellow"
  gsn_text(wks,plot,"DMS",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="red"
  gsn_text(wks,plot,"SO2",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="purple"
  gsn_text(wks,plot,"SO3",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="blue"
  gsn_text(wks,plot,"H2SO4",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="orange"
  gsn_text(wks,plot,"SO4",xpos,ypos,txres)

  draw(plot)
  frame(wks)
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
    year=1993
  end if
  print(year)

  plotS(rootpath, caseid, year)
end if
