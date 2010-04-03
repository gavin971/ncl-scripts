; Program checks solar lean netCDF file that integral of spectrum
; equals TSI, and that TSI is same as for file
; SOLAR_TSI_Lean_1610-2140_annual_c100301.nc

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

ncf = addfile("spectra_1610_2009a_31Mar10.nc","r")

ssi=ncf->ssi
tsi=ncf->tsi
band_width=ncf->band_width
time=ncf->time
date=ncf->date

tsi_sum = ssi # band_width
tsi_sum = tsi_sum/1000.d0 ; Convert mW to W
tsi_sum!0="time"
tsi_sum&time=time
tsi_diff=tsi-tsi_sum

tsi_sum2 = ssi(:,1:) # band_width(1:)
tsi_sum2 = tsi_sum2/1000.d0 ; Convert mW to W
tsi_sum2!0="time"
tsi_sum2&time=time
tsi_diff3=tsi-tsi_sum2

ncf2 = addfile("SOLAR_TSI_Lean_1610-2140_annual_c100301.nc","r")
tsi2= ncf2->tsi

tsi_diff2=tsi-tsi2(0:dimsizes(tsi)-1)


ncf3 = addfile("spectra_1610_2000a_21Jan09.nc","r")
tsi3=ncf3->tsi

year=round(date/10000.d0,0)+0.5d0

do i = 220,260
  print(sprintf("%6.1f",year(i))+"  "+sprintf("%16.11f",tsi(i))+\
        "  "+sprintf("%16.11f",tsi2(i))+"  "+sprintf("%16.11f",tsi_diff2(i)))
end do

do i = 389,390
  print(sprintf("%6.1f",year(i))+"  "+sprintf("%16.11f",tsi(i))+\
        "  "+sprintf("%16.11f",tsi2(i))+"  "+sprintf("%16.11f",tsi_diff2(i))+\
        "  "+sprintf("%16.11f",tsi3(i)))
end do
do i = 391,399
  print(sprintf("%6.1f",year(i))+"  "+sprintf("%16.11f",tsi(i))+\
        "  "+sprintf("%16.11f",tsi2(i))+"  "+sprintf("%16.11f",tsi_diff2(i)))
end do


wks = gsn_open_wks("pdf", "tsi_diff")

  res                   = True                      ; plot mods desired
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@vpHeightF             = 0.5              ; change aspect ratio
  res@vpWidthF              = 1.0

  res@trXMinF = 1610.d0
  res@trXMaxF = 2010.d0
  res@tmXBLabelFontHeightF  = 0.02
;  res@tiXAxisString="Time (years)"
  res@tiYAxisString="TSI - SSI sum (W/m~S~2~N~)"
  
  plot = new(2, graphic)
  data=new((/2,dimsizes(tsi)/),double)
  data(0,:)=tsi_diff
  data(1,:)=tsi_diff3
  plot(0) = gsn_xy(wks, year, tsi_diff, res)
  
  res@tiYAxisString="TSI new - TSI old (W/m~S~2~N~)"
  plot(1) = gsn_xy(wks, year, tsi_diff2, res)

  resP                      = True
  resP@gsnMaximize          = True
;  resP@gsnPanelLabelBar     = True
;  resP@lbTitleOn            =  True                ; turn on title
;  resP@lbTitlePosition      = "Right"              ; title position
;  resP@lbTitleFontHeightF   = .015                ; make title smaller
;  resP@lbTitleDirection     = "Across"             ; title direction
  resP@gsnPanelXWhiteSpacePercent = 5
  gsn_panel(wks, plot, (/ 2, 1 /), resP)
  
delete(wks)
delete(res)
  
wks = gsn_open_wks("pdf", "tsi_vs_ssi_1830-1870")
  res                   = True                      ; plot mods desired
;  res@gsnDraw           = False
;  res@gsnFrame          = False
;  res@vpHeightF             = 0.5              ; change aspect ratio
;  res@vpWidthF              = 1.0
;  res@trXMinF = 1610.d0
;  res@trXMaxF = 2010.d0


  res@trXMinF = 1830.d0
  res@trXMaxF = 1870.d0
  res@trYMinF = 1360.4
  res@trYMaxF = 1361.6
  res@tmXBLabelFontHeightF  = 0.015
  res@gsnMaximize          = True
  res@gsnPaperOrientation = "landscape"
  res@tmXBMinorPerMajor = 10

;  res@tmYBLabelFontHeightF  = 0.01
  res@tiYAxisFontThicknessF = 0.1
;  res@tiXAxisString="Time (years)"
  res@tiYAxisString="TSI (W/m~S~2~N~)"
  res@xyLineColors      = (/"blue","red","green"/)
  res@xyDashPatterns      = (/0,0,1/)
  
  delete(data)
  data=new((/2,dimsizes(tsi)/),double)
  data(0,:)=tsi
  data(1,:)=tsi_sum
;  data(2,:)=tsi_sum2
  plot1 = gsn_xy(wks, year, data, res)
  
;  delete(data)
;  delete(res@xyLineColors)
;  res@xyLineColors      = (/"red","green"/)
;  delete(res@xyDashPatterns)
;  res@xyDashPatterns      = (/0,1/)
;  res@tiYAxisString="TSI - SSI sum (W/m~S~2~N~)"
;  data=new((/2,dimsizes(tsi)/),double)
;  data(0,:)=tsi_diff
;  data(1,:)=tsi_diff3
;  plot(1) = gsn_xy(wks, year, tsi_diff, res)



; tsi_rat=tsi/tsi_sum
;  delete(res@xyLineColors)
;  res@xyLineColors      = (/"black"/)
;  res@tiYAxisString="TSI / SSI sum"
  
;  plot(1) = gsn_xy(wks, year, tsi_rat, res)
  
  
;  resP                      = True
;  resP@gsnMaximize          = True


;  resP@gsnPaperOrientation = "landscape"
;  resP@gsnPanelLabelBar     = True
;  resP@lbTitleOn            =  True                ; turn on title
;  resP@lbTitlePosition      = "Right"              ; title position
;  resP@lbTitleFontHeightF   = .015                ; make title smaller
;  resP@lbTitleDirection     = "Across"             ; title direction


;  resP@gsnPanelXWhiteSpacePercent = 5
;  gsn_panel(wks, plot, (/ 2, 1 /), resP)
 
