; Ice Mass - total and by Latitude and Level

;load "tests/common.ncl"
;load "tests/gsn_csm_pres_hgt_fix.ncl" 

undef("plotS")
procedure plotS(suffix, WACCM_RUN, startFile, stopFile)

begin

  ; Open the histories
  files = systemfunc("ls " + WACCM_RUN + "/atm/hist/*h0*.nc")
  
  history = addfiles(files(startFile:stopFile), "r")
  
  z     = history[:]->Z3
  z     = (/ z / 1000. /)

  ; Total the Gases
  h2so4 = history[:]->H2SO4
  hso3  = history[:]->HSO3
  ocs   = history[:]->OCS
  s     = history[:]->S
  so    = history[:]->SO
  so2   = history[:]->SO2
  so3   = history[:]->SO3
  
  totalgas = h2so4
  totalgas = (/ totalgas + hso3 + ocs + s + so + so2 + so3 /)
  
  
  ; Add in the pure sulfates
  pure  = history[:]->PURSULMR
  pure  = (/ pure * 28.89 / 98. /)
  
  ; For the mixed sulfates, the core mass needs to be subtracted.
  mix   = history[:]->MIXSULMR
  
  nBins = 28
  do i = 1, nBins
    fieldname = "MIXCOR" + sprinti("%2.2i", i)
    core = history[:]->$fieldname$
    
    mix = (/ mix - core /)
  end do
  
  mix = (/ mix * 28.89 / 98. /)
  mix = (/ where(mix .le. 0., 0., mix) /)
  
  totalaer = pure
  totalaer = (/ pure + mix /)
  
  ; Calculate the total in pptv
  totals = totalgas
  totals = (/ totalgas + totalaer /)
  totals = (/ totals * 1e12 /)
  totals@units = "pptv"
  totals@long_name = "Total Sulfur"

  ; Make a zonal and time average.
  tmp = totals(lev|:, lat|:, time|:, lon|:)
  ts  = dim_avg_Wrap(dim_avg_Wrap(tmp))
  gts = dim_avg_Wrap(ts)
  delete(tmp)

  tmp = z(lev|:, lat|:, time|:, lon|:)
  gz  = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(tmp)))
  delete(tmp)

  tmp = ocs(lev|:, lat|:, time|:, lon|:)
  tmp = (/ tmp * 1e12 /)
  gocs = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(tmp)))
  delete(tmp)
   
  tmp = so2(lev|:, lat|:, time|:, lon|:)
  tmp = (/ tmp * 1e12 /)
  gso2 = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(tmp)))
  delete(tmp)
   
  tmp = so3(lev|:, lat|:, time|:, lon|:)
  tmp = (/ tmp * 1e12 /)
  gso3 = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(tmp)))
  delete(tmp)
   
  tmp = h2so4(lev|:, lat|:, time|:, lon|:)
  tmp = (/ tmp * 1e12 /)
  gh2so4 = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(tmp)))
  delete(tmp)
   
  tmp = pure(lev|:, lat|:, time|:, lon|:)
  tmp = (/ tmp * 1e12 /)
  gpure = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(tmp)))
  delete(tmp)
   
  tmp = mix(lev|:, lat|:, time|:, lon|:)
  tmp = (/ tmp * 1e12 /)
  gmix = dim_avg_Wrap(dim_avg_Wrap(dim_avg_Wrap(tmp)))
  delete(tmp)
  
 
  ; Open a workstation
  wks = gsn_open_wks("pdf", "images/sulfur_vmr" + suffix)
;  gsn_define_colormap(wks,"BlWhRe")       
      
  plot = new(2, graphic)

  resP                       = True
  resP@gsnMaximize           = True
;  resP@gsnPaperOrientation   = "portrait"
  resP@gsnPanelLabelBar     = False

;  resP@lbTitleOn            = True                ; turn on title
;  resP@lbTitlePosition      = "Right"              ; title position
;  resP@lbTitleFontHeightF   = .015                ; make title smaller
;  resP@lbTitleDirection     = "Across"             ; title direction
;  resP@lbTitleString = "ppmv"
;  resP@lbLabelStride        = 2
;  resP@txString = "Water Vapor Anomaly, 20N-20S"
;  resP@gsnPanelXWhiteSpacePercent = 5.
;  resP@gsnPanelYWhiteSpacePercent = 5.

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

  ts&lev@units = "hPa"
  ts&lev@long_name = "Pressure (hPa)"

print(min(ts({minP:maxP},:)))
print(max(ts({minP:maxP},:)))
  plot(0) = gsn_csm_pres_hgt(wks, ts, res)
  

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

  plot(1) = gsn_csm_xy(wks, gts, gz, resxy)

  resxy@xyLineThicknesses     = (/ 4.0 /)
  resxy@xyLineColors          = (/ "green" /)
  overlay(plot(1), gsn_csm_xy(wks, gocs, gz, resxy))

  resxy@xyLineColors          = (/ "red" /)
  overlay(plot(1), gsn_csm_xy(wks, gso2, gz, resxy))

  resxy@xyLineColors          = (/ "purple" /)
  overlay(plot(1), gsn_csm_xy(wks, gso3, gz, resxy))

  resxy@xyLineColors          = (/ "blue" /)
  overlay(plot(1), gsn_csm_xy(wks, gh2so4, gz, resxy))

  resxy@xyLineColors          = (/ "green" /)
  overlay(plot(1), gsn_csm_xy(wks, gocs, gz, resxy))

  resxy@xyLineColors          = (/ "orange" /)
  overlay(plot(1), gsn_csm_xy(wks, gmix + gpure, gz, resxy))
  resxy@xyDashPatterns        = (/ 2 /)
  overlay(plot(1), gsn_csm_xy(wks, gmix, gz, resxy))
  resxy@xyDashPatterns        = (/ 1 /)
  overlay(plot(1), gsn_csm_xy(wks, gpure, gz, resxy))

  gsn_panel(wks, plot, (/  1, 2 /), resP)
end

WACCM_RUN = getenv("WACCM_RUN")

if (ismissing(WACCM_RUN)) then
  WACCM_RUN = "results"
end if
print("Using WACCM_RUN: " + WACCM_RUN)
print("")

plotS("", "SD2004_pmc_81grf_4jso31so20nn", 12, 23)
plotS("_n", "SD2004_pmc_81grf_4jso31so20", 12, 23)
plotS("_n.001", "SD2004_pmc_81grf_4jso31so20.001", 12, 23)
plotS("_n.1", "SD2004_pmc_81grf_4jso31so20.1", 12, 23)
