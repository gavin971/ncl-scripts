; Meteoritic Dust
;
; Monthly average plots

; Constants
;
; NOTE: These are from CAM, shr_const_mod.F90.
PI      = 3.14159265358979323846  ; pi
G       = 9.80616                 ; acceleration of gravity ~ m/s^2
BOLTZ   = 1.38065e-23             ; Boltzmann's constant ~ J/K/molecule
AVOGAD  = 6.02214e26              ; Avogadro's number ~ molecules/kmole
RGAS    = AVOGAD*BOLTZ            ; Universal gas constant ~ J/K/kmole
MWDAIR  = 28.966                  ; molecular weight dry air ~ kg/kmole
RDAIR   = RGAS/MWDAIR             ; Dry air gas constant ~ J/K/kg
MWH2O   = 18.000                  ; molecular weight of water vapor
PTOP    = 100.0                   ; pressure at the top of the model (non-zero) ~ Pa
REARTH  = 6.37122e6               ; Earth's radius (m)

load "ncl-scripts/core.ncl"

;  Set up the particle bins.
;
;  For each particle group, the mass of a particle in
;  bin i is <rmrat> times that in bin i-1
;
;    rmass(NBIN)     =  bin center mass [g]
;    r(NBIN)         =  bin mean (volume-weighted) radius [cm]
;    vol(NBIN)       =  bin center volume [cm^3]
;    dr(NBIN)        =  bin width in radius space [cm]
;    dv(NBIN)        =  bin width in volume space [cm^3]
;    dm(NBIN)        =  bin width in mass space [g]
;
; NOTE: This logic is from CARMA setupbins.f.

undef("setupbins")
procedure setupbins(nbin:integer, rmin:float, rmrat:float, rho:float, rmass[*]:float, r[*]:float, dr[*]:float)

local cpi, i, vrfact

begin

  ; Allocate arrays for the return variables.
  rmass = new(nbin, float)
  r = new(nbin, float)
  dr = new(nbin, float)
  
  r@units = "cm"
  r@long_name = "Radius"

  ; NOTE: These aren't currently being returned, but are calculated
  ; in setupbins.
  rmassup = new(nbin, float)
  dm = new(nbin, float)
  vol = new(nbin, float)
  rup = new(nbin, float)
  rlow = new(nbin, float)
  
  ; The first bin is based upon rmin.
  cpi = 4./3. * PI
  rmass(0) = cpi * rho * (rmin^3)
  vrfact = ((3./2. / PI / (rmrat + 1.))^(1./3.)) * (rmrat^(1./3.) - 1.)

  ; The remaining bins have a mass that is rmrat larger than the
  ; previous bin.
  do i = 0, nbin-1

    rmass(i)   = rmass(0) * (rmrat^(i))
    rmassup(i) = 2. * rmrat / (rmrat + 1.) * rmass(i)
    dm(i)      = 2. * (rmrat - 1.) / (rmrat + 1.) * rmass(i)
    vol(i)     = rmass(i) / rho
    r(i)       = (rmass(i) / rho / cpi) ^ (1./3.)
    rup(i)     = (rmassup(i) / rho / cpi) ^ (1./3.)
    dr(i)      = vrfact * ((rmass(i)/rho) ^ (1./3.))
    rlow(i)    = rup(i) - dr(i)
  end do
end


; Calculate parameters of the air based upon the data in the history file.
;
undef("setupair")
procedure setupair(historyFile:file, P, arho)

local hyam, hybm, P0, PS, Q, T, arho

begin
  hyam 	= historyFile->hyam
  hybm 	= historyFile->hybm
  P0 	= historyFile->P0
  PS 	= historyFile->PS
  Q     = historyFile->Q
  T     = historyFile->T
 
  ; Calculate the air pressures at the interfaces.
  P = pres_hybrid_ccm(PS(:,:,:), P0, hyam(:) ,hybm(:))
  
  ; Calculate the density assuming an ideal gas.
  arho(:,:,:,:)  = P(:,:,:,:) / (RDAIR * T(:,:,:,:))
end

; Get the number density for each bin.
undef("getN")
procedure getN(history:file, shortName:string, rmass:float, arho:float, N:float)
  
local b, binName, nBins, mmr

begin
  
  ; Iterate over the bins
  nBins = dimsizes(rmass)
  
  do b = 0, nBins-1
    binName = shortName + sprinti("%2.2i", b+1)
    
    ; Get the mass mixing ratio.
    mmr = history->$binName$
    
    ; Compute the particle concentration.
    N(:, b, :, :, :) = (mmr(lev|:,time|:,lat|:,lon|:) * arho(lev|:,time|:,lat|:,lon|:)) / rmass(b)    
  end do
end 


undef("monave_base")
procedure monave_base(WACCM_RUN:string)

local files, nFiles, historyFile

begin

  ; Read in the global average data for one year.
  history = addfile("analysis/" + WACCM_RUN + "/h0.2034.mon.nc", "r")

  ; Compute the pressure levels ...
  lev     = history->lev
  nLevs   = dimsizes(lev)
  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)"
  
  ; Determine an approximate dz/dp.
  hyam 	= history->hyam
  hybm 	= history->hybm
  P0 	= history->P0
  PS 	= history->PS
  Q     = history->Q
  T     = history->T
 
  ; Calculate the air pressures at the interfaces.
  P = pres_hybrid_ccm(PS(:,:,:), P0, hyam(:) ,hybm(:))
  
  ; Calculate the density assuming an ideal gas.
  arho  = P(:,:,:,:) / (RDAIR * T(:,:,:,:))
  dzdp = - 1 / (arho * G)
  
  
  ; Make a 12 plot panel ...
  nPlots = 12
  plot = new((/ nPlots /), graphic)
  
  resP                      = True
  resP@gsnMaximize          = True
;  resP@gsnPaperOrientation  = "portrait"
  resP@gsnPanelLabelBar     = True
        
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@lbLabelBarOn          = False
  res@tmYLMode   = "Explicit" ; Define own tick mark labels.
  res@tmYLValues = (/1000.,  850., 700., 500., 400., 300., 250.,\
                         200., 150., 100.,  70., 50., 30., 10., 1., 0.1, 0.01, 0.001, 0.0001/)
  res@tmYLLabels = (/"1000","850","700","500","400","300","250",\
                         "200","150","100", "70", "50", "30", "10", "1", "0.1", "0.01", "0.001", "0.0001"/)
  res@tmYLMinorOn= False        ; No minor tick marks.
  res@cnFillOn              = True                  ; turn on color fill
  res@cnInfoLabelOn         = False
  res@gsnSpreadColors       = True                  ; use full range of colors

  months = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" /)

  ; Open a workstation
  wks = gsn_open_wks("pdf","images/" + WACCM_RUN + "/monave_base")
  gsn_define_colormap(wks,"BlWhRe")                ; choose colormap
  
  zonal = dim_avg_Wrap(history->U)
  zonal&lev = lev_p
  
  resP@txString             = "Zonal Wind"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 10.0                   ; contour interval
  res@cnMinLevelValF       = -100.0             ; min level
  res@cnMaxLevelValF       = 100.0            ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)


  zonal = dim_avg_Wrap(history->V)
  zonal&lev = lev_p
  
  resP@txString             = "Meridional Wind"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 5.0                   ; contour interval
  res@cnMinLevelValF       = -40.0             ; min level
  res@cnMaxLevelValF       = 40.0            ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)


  zonal = dim_avg_Wrap(history->OMEGA * dzdp) * 100.0
  zonal&lev = lev_p
  zonal@long_name = "Vertical Wind"
  zonal@units = "cm/s"
  
  
  resP@txString             = "Vertical Wind"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 0.5           ; contour interval
  res@cnMinLevelValF       = -5.0           ; min level
  res@cnMaxLevelValF       = 5.0           ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)


  gsn_define_colormap(wks,"rainbow")                ; choose colormap
  
  zonal = dim_avg_Wrap(history->T)
  zonal&lev = lev_p
  
  resP@txString             = "Temperature"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 20.0                   ; contour interval
  res@cnMinLevelValF       = 120.0             ; min level
  res@cnMaxLevelValF       = 300.0            ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)


  zonal = dim_avg_Wrap(history->Q) * 1e6 * (MWDAIR / MWH2O)
  zonal&lev = lev_p
  zonal@long_name = "Specific Humidity"
  zonal@units = "ppmv"
  
  resP@txString             = "Specific Humidity"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = .5                  ; contour interval
  res@cnMinLevelValF       = 0            ; min level
  res@cnMaxLevelValF       = 7             ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)
end


undef("monave_dust")
procedure monave_dust(WACCM_RUN:string)

local files, nFiles, historyFile

begin

  ; Read in the global average data for one year.
  history = addfile("analysis/" + WACCM_RUN + "/h0.2034.mon.nc", "r")

  ; Compute the pressure levels ...
  lev     = history->lev
  nLevs   = dimsizes(lev)
  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)"
  
  
  
  ; Make a 12 plot panel ...
  nPlots = 12
  plot = new((/ nPlots /), graphic)
  
  resP                      = True
  resP@gsnMaximize          = True
;  resP@gsnPaperOrientation  = "portrait"
  resP@gsnPanelLabelBar     = True
        
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@lbLabelBarOn          = False
  res@tmYLMode   = "Explicit" ; Define own tick mark labels.
  res@tmYLValues = (/1000.,  850., 700., 500., 400., 300., 250.,\
                         200., 150., 100.,  70., 50., 30., 10., 1., 0.1, 0.01, 0.001, 0.0001/)
  res@tmYLLabels = (/"1000","850","700","500","400","300","250",\
                         "200","150","100", "70", "50", "30", "10", "1", "0.1", "0.01", "0.001", "0.0001"/)
  res@tmYLMinorOn= False        ; No minor tick marks.
  res@cnFillOn              = True                  ; turn on color fill
  res@cnInfoLabelOn         = False
  res@gsnSpreadColors       = True                  ; use full range of colors

  months = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" /)

  ; Open a workstation
  wks = gsn_open_wks("pdf","images/" + WACCM_RUN + "/monave_dust")
  gsn_define_colormap(wks,"rainbow")                ; choose colormap
  
  zonal = dim_avg_Wrap(history->DUSTND)
  zonal&lev = lev_p
  zonal = dim_avg_Wrap(history->DUSTND) * 1e-3
  zonal@units = "10~S~3~N~ /cm~S~3~N~"
  zonal@long_name = "Dust Concentration"
  
  resP@txString             = "Dust Concentration"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 0.5                  ; contour interval
  res@cnMinLevelValF       = 0.0             ; min level
  res@cnMaxLevelValF       = 7.0            ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)


  zonal = dim_avg_Wrap(history->DUSTAD)
  zonal&lev = lev_p
  
  resP@txString             = "Dust Surface Area Density"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = .1                   ; contour interval
  res@cnMinLevelValF       = 0             ; min level
  res@cnMaxLevelValF       = 1.6            ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)


  zonal = dim_avg_Wrap(history->DUSTMD) * 1e15
  zonal&lev = lev_p
  zonal@units = "10~S~-15~N~ g/cm~S~3~N~"
  zonal@long_name = "Dust Density"
  
  
  resP@txString             = "Dust Density"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 5.0           ; contour interval
  res@cnMinLevelValF       = 0.0          ; min level
  res@cnMaxLevelValF       = 60.0           ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)


  zonal = dim_avg_Wrap(history->DUSTRE) * 1000.0
  zonal&lev = lev_p
  zonal@units = "nm"
  zonal@long_name = "Dust Effective Radius"
  
  resP@txString             = "Dust Effective Radius"
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 5.0                   ; contour interval
  res@cnMinLevelValF       = 0.0             ; min level
  res@cnMaxLevelValF       = 70.0            ; max level
  print(min(zonal(:, {0.0:150.0}, :)))
  print(max(zonal(:, {0.0:150.0}, :)))
  
  do i = 0, nPlots-1
    res@tiMainString          = months(i)
    plot(i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:150.0}, :), res)
  end do
  
  gsn_panel(wks, plot, (/ 4, 3 /), resP)
end



WACCM_RUN = getenv("WACCM_RUN")

if (ismissing(WACCM_RUN)) then
  WACCM_RUN = "results"
end if
print("Using WACCM_RUN: " + WACCM_RUN)
print("")

monave_base(WACCM_RUN)
monave_dust(WACCM_RUN)
