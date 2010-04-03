; Program plots a contour of solar spectral irradiance vs.
; time (years 1610-2000) and wavelength bins

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

ncf = addfile("spectra_1610_2009a_31Mar10.nc","r")

ssi=ncf->ssi
date=ncf->date
year=round(date/10000.d0,0)+0.5d0
wvl=ncf->wvl

wks = gsn_open_wks("pdf", "contour_ssi")

data = log(ssi)
data!0=ssi!0
data!1=ssi!1
data&time = year
data&wavelength=wvl
;data&wavelength=log10(wvl)
print("max: "+max(data))
print("min: "+min(data))

gsn_define_colormap(wks,"gui_default")                ; choose colormap

res = True
res@gsnMaximize            = True
res@cnFillOn              = True                  ; turn on color fill
res@gsnSpreadColors       = True                  ; use full range of colors
;res@trYLog            = True                       ; log Y-axis
;res@trYMinF            = 120.d0
;res@trYMaxF            = 100000.d0

  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 1.0                   ; contour interval
  res@cnMinLevelValF       = -10.0             ; min level
  res@cnMaxLevelValF       = 8.0            ; max level
  res@cnLineThicknessF     = 2.0                   ; color of second contours
  res@lbLabelAngleF       = 90

plot = gsn_csm_contour(wks, data(wavelength|:,time|:), res)

  res@cnLevelSpacingF      = 1.0                   ; contour interval
  res@cnMinLevelValF       = -4.0             ; min level
  res@cnMaxLevelValF       = 8.0            ; max level
plot = gsn_csm_contour(wks, data(wavelength|0:500,time|:), res)

  res@cnLevelSpacingF      = 0.5                   ; contour interval
  res@cnMinLevelValF       = -4.0             ; min level
  res@cnMaxLevelValF       = 4.0            ; max level
plot = gsn_csm_contour(wks, data(wavelength|0:100,time|:), res)

  res@cnLevelSpacingF      = 0.5                   ; contour interval
  res@cnMinLevelValF       = -4.0             ; min level
  res@cnMaxLevelValF       = 2.5            ; max level
plot = gsn_csm_contour(wks, data(wavelength|0:10,time|:), res)
