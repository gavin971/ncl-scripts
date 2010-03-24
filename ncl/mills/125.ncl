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
CPDAIR  = 1.00464e3               ; specific heat of dry air ~ J/kg/K
MWH2O   = 18.000                  ; molecular weight of water vapor
PTOP    = 100.0                   ; pressure at the top of the model (non-zero) ~ Pa
REARTH  = 6.37122e6               ; Earth's radius (m)

load "ncl-scripts/common.ncl"

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


undef("monave_za")
procedure monave_za()

local files, nFiles, historyFile

begin
  print("")
  print("  Plot monthly average dust densities ...")
  
  ; Set up the particle bins ...
  nBins = 28
  rmin  = 2e-8 ; cm
  rmrat = 2.0
  rho   = 2.0  ; g/cm-3

  rmass = new(nBins, float)
  r = new(nBins, float)
  dr = new(nBins, float)

  print("")
  print("  Assuming ...")
  print("    nBins    =  " + nBins)
  print("    rmin     =  " + rmin)
  print("    rmrat    =  " + rmrat)
  print("    rho      =  " + rho)
  
  setupbins(nBins, rmin, rmrat, rho, rmass, r, dr)
  
  ; Open the control history file
  history = addfile("/home/bardeen/waccm/carma/analysis/125_lt_gw2_28/h0.2028.mon.nc", "r")
  
  ; Compute the pressure levels ...
  time     = history->time
  nTimes   = dimsizes(time)
  lev     = history->lev
  nLevs   = dimsizes(lev)
  lat     = history->lat
  nLats   = dimsizes(lat)
  lon     = history->lon
  nLons   = dimsizes(lon)
  
  ; Determine an approximate dz/dp.
  hyam 	= history->hyam
  hybm 	= history->hybm
  hyai 	= history->hyai
  hybi 	= history->hybi
  P0 	  = history->P0
  PS 	  = history->PS
  Q     = history->Q
  T     = history->T
  
  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)"
  
  ; Make three sets of data:
  ;  control
  ;  test
  ;  difference
  month = (/ 1, 7 /)
  nMonths = dimsizes(month)
  runs = (/ "January", "July" /)
  nRuns = dimsizes(runs)
  
  Ts = new((/ nRuns, nLevs, nLats /), float)
  Ts!0 = "run"
  Ts&run = runs
  Ts!1 = "lev"
  Ts&lev = lev_p
  Ts!2 = "lat"
  Ts&lat = lat
  Ts@long_name = "Temperature"
  Ts@units = "K"
  
  N1s = new((/ nRuns, nLevs, nLats /), float)
  N1s!0 = "run"
  N1s&run = runs
  N1s!1 = "lev"
  N1s&lev = lev_p
  N1s!2 = "lat"
  N1s&lat = lat
  N1s@long_name = "Concentration, r >= 1 nm"
  N1s@units = "# cm~S~-3~N~"  

  NDs = new((/ nRuns, nLevs, nLats /), float)
  NDs!0 = "run"
  NDs&run = runs
  NDs!1 = "lev"
  NDs&lev = lev_p
  NDs!2 = "lat"
  NDs&lat = lat
  NDs@long_name = "Concentration"
  NDs@units = "# cm~S~-3~N~"  

  ADs = new((/ nRuns, nLevs, nLats /), float)
  ADs!0 = "run"
  ADs&run = runs
  ADs!1 = "lev"
  ADs&lev = lev_p
  ADs!2 = "lat"
  ADs&lat = lat
  ADs@long_name = "Surface Area Density"
  ADs@units = "um~S~2~N~ cm~S~-3~N~"  

  MDs = new((/ nRuns, nLevs, nLats /), float)
  MDs!0 = "run"
  MDs&run = runs
  MDs!1 = "lev"
  MDs&lev = lev_p
  MDs!2 = "lat"
  MDs&lat = lat
  MDs@long_name = "Mass Density"
  MDs@units = "x10~S~-16~N~ g cm~S~-3~N~"  

  REs = new((/ nRuns, nLevs, nLats /), float)
  REs!0 = "run"
  REs&run = runs
  REs!1 = "lev"
  REs&lev = lev_p
  REs!2 = "lat"
  REs&lat = lat
  REs@long_name = "Effective Radius"
  REs@units = "nm"  

  
  ; Calculate the air pressures at the interfaces.
  print("")
  print("  Calculating air pressure and density ...")
  P = Q
  P = (/ pres_hybrid_ccm(PS(:,:,:), P0, hyam(:) ,hybm(:)) /)
  P@units         = "Pa"
  P@long_name     = "pressure"
  
  ; Calculate the density assuming an ideal gas.
  arho            = P
  arho@units      = "g/cm3"
  arho@long_name  = "air density"

  setupair(history, P, arho)
  arho = arho(:,:,:,:) * 1e3 / 1e6
  
  ; Get the Number density (#/cm-3) for each bin at the 
  ; specified time and location.
  print("")
  print("  Calculating dust number, area and mass densities ...")
  N = new((/ nLevs, nBins, nTimes, nLats, nLons /), float)
  N!0 = "lev"
  N&lev = lev
  N!1 = "radius"
  N&radius = r
  N!2 = "time"
  N&time = time
  N!3 = "lat"
  N&lat = lat
  N!4 = "lon"
  N&lon = lon
  N@units         = "#/cm3"
  N@long_name     = "number density"
  getN(history, "DUST", rmass, arho, N)


  ; Calculate the zonal average of total number, area and mass in the
  ; desired bin range.
  firstBin  =  7           ; smallest particle to include in the totals
  lastBin   = nBins-1     ; largest particle to include in the totals
  print("  smallest radius = " + r(firstBin))

  N1 = dim_sum_Wrap(N(time|:, lev|:, lat|:, lon|:, radius|firstBin:lastBin))
  ND = dim_sum_Wrap(N(time|:, lev|:, lat|:, lon|:, radius|:))
  AD = ND
  MD = ND
  RE = ND

  AD = (/ 0.0 /)
  MD = (/ 0.0 /)
  RE = (/ RE@_FillValue /)

  do b = 0, nBins-1
    NS = N(time|:, lev|:, lat|:, lon|:, radius|b)
    AD(:, :, :, :) = AD(:, :, :, :) + (NS * 4.0 * PI * (r(b)^2) * 1e8)
    MD(:, :, :, :) = MD(:, :, :, :) + (NS * rmass(b))
  end do

  ; For the effective radius calculaton, we need to do a divsion that
  ; might have zero values in it.
  ND1 = ndtooned(ND)
  MD1 = ndtooned(MD)
  AD1 = ndtooned(AD)
  RE1 = ndtooned(RE)

  iad = ind(ND1.gt.1.0)

  RE1(iad) = (/ 3.0 * (MD1(iad) / (AD1(iad) * rho / 1e8)) * 1e4 /)

  RE = (/ onedtond(RE1, dimsizes(RE)) /)

  delete(RE1)
  delete(AD1)
  delete(MD1)
  delete(ND1)
  delete(iad)

  do i = 0, nMonths - 1
    Ts(i, :, :) = dim_avg(T(month(i), :, :, :))
    N1s(i, :, :) = dim_avg(N1(month(i), :, :, :))
    NDs(i, :, :) = dim_avg(ND(month(i), :, :, :))
    ADs(i, :, :) = dim_avg(AD(month(i), :, :, :))
    MDs(i, :, :) = dim_avg(MD(month(i), :, :, :)) * 1e16
    REs(i, :, :) = dim_avg(RE(month(i), :, :, :)) * 1000.0
  end do
      
  delete(N1)
  delete(ND)
  delete(AD)
  delete(MD)
  delete(RE)
  delete(N)

  ; Make a 4 plot panel ...
  print("")
  print("  Plotting ...")
 
  minP      = 0.0001        ; largest pressure to plot ...
  maxP      = 0.1        ; largest pressure to plot ...
  
  nPlots = 2
  plot = new((/ nPlots /), graphic)
  plotT = new((/ nPlots /), graphic)
  
  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnPanelLabelBar     = True
;  resP@lbTitleOn            =  True                ; turn on title
;  resP@lbTitlePosition      = "Right"              ; title position
;  resP@lbTitleFontHeightF   = .015                ; make title smaller
;  resP@lbTitleDirection     = "Across"             ; title direction
  resP@gsnPanelXWhiteSpacePercent = 0.5
  resP@gsnPanelYWhiteSpacePercent = 0.5
  resP@lbLabelStride        = 1
        
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@lbLabelBarOn          = False
  res@cnLinesOn             = False
  res@cnLinesOn             = False
  res@cnFillOn              = True                  ; turn on color fill
  res@cnInfoLabelOn         = False
  res@cnLineLabelsOn        = False
  res@gsnSpreadColors       = True                  ; use full range of colors
  res@cnFixFillBleed        = True
;  res@cnFillMode            = "CellFill"
;  res@lbTitlePosition       = "Right"              ; title position
;  res@lbTitleDirection      = "Across"             ; title direction
;  resP@lbTitleFontHeightF   = .015                ; make title smaller
  res@tmXBLabelFontHeightF  = 0.03   ; increase tickmark font size
  res@tmYRLabelFontHeightF  = 0.03   ; increase tickmark font size
  res@tmYLLabelFontHeightF  = 0.03   ; increase tickmark font size
  res@tmYLLabelAngleF      = 60
  
  resO                      = True
  resO@gsnDraw              = False
  resO@gsnFrame             = False
  resO@lbLabelBarOn         = False
  resO@cnLineColor          = "Black"                   ; color of second contours
  resO@gsnRightString       = ""    
  resO@gsnLeftString        = ""
  resO@gsnCenterString      = ""
  resO@cnLevelSelectionMode = "ManualLevels" 
  resO@cnLevelSpacingF      = 10.0                   ; contour interval
  resO@cnMinLevelValF       = 110.0             ; min level
  resO@cnMaxLevelValF       = 160.0            ; max level
  resO@cnLineThicknessF     = 2.0                   ; color of second contours

  months = (/ "Dec", "January", "Feb", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov" /)

  ; Open a workstation
  wks = gsn_open_wks("pdf","images/125")
  gsn_define_colormap(wks,"gui_default")                ; choose colormap
  
  print("")
  print("    temperature overlay")
  zonal = Ts
  
  do i = 0, nPlots-1
    plotT(i) = gsn_csm_contour(wks, zonal(i, {0.0:maxP}, :), resO)          ; draw second plot
  end do

  print("")
  print("    number density")
  zonal                     = N1s
  
;  res@lbTitleString        = zonal@units
;  res@lbLabelStride         = 3
  res@cnLevelSelectionMode  = "ManualLevels" 

  print("      min :  " + min(zonal(:, {minP:maxP}, :)))
  print("      max :  " + max(zonal(:, {minP:maxP}, :)))
  res@cnLevelSpacingF       = 5.0                  ; contour interval
  res@cnMinLevelValF        = 5.0                  ; contour interval
  res@cnMaxLevelValF        = 100.0                 ; contour interval  
  res@tiMainString          = runs(0)
  res@tiMainFontHeightF         = 0.04
  
  res@tiMainString          = months(month(0))
  plot(0) = gsn_csm_pres_hgt(wks, zonal(0, {minP:maxP}, {-90.0:-30.0}), res)
  overlay(plot(0), plotT(0))

  res@tiMainString          = months(month(1))
  plot(1) = gsn_csm_pres_hgt(wks, zonal(1, {minP:maxP}, {30.0:90.0}), res)
  overlay(plot(1), plotT(1))
    

;  print("")
;  print("    mass density")
;  zonal = MDs
;  zonal&lev = lev_p

;  res@lbTitleString        = zonal@units
;  print("      min :  " + min(zonal(:, {0.0:maxP}, :)))
;  print("      max :  " + max(zonal(:, {0.0:maxP}, :)))

;  res@cnLevelSelectionMode = "ExplicitLevels"
;  res@cnLevels = (/ .2, .4, .8, 1, 2, 4, 8, 10, 20, 40, 80, 100, 200 /)
;  res@lbLabelStride        = 1

;  res@tiMainString          = ""
;  do i = 0, nRuns-1
;    plot(2+i) = gsn_csm_pres_hgt(wks, zonal(i, {0.0:maxP}, :), res)
;  end do

  gsn_panel(wks, plot, (/ 1, 2 /), resP)

end



undef("monave_polar")
procedure monave_polar()

local files, nFiles, historyFile

begin
  print("")
  print("  Plot monthly average dust densities ...")
  
  print("")
  print("  Plot monthly average dust densities ...")
  
  ; Set up the particle bins ...
  nBins = 28
  rmin  = 2e-8 ; cm
  rmrat = 2.0
  rho   = 2.0  ; g/cm-3

  rmass = new(nBins, float)
  r = new(nBins, float)
  dr = new(nBins, float)

  print("")
  print("  Assuming ...")
  print("    nBins    =  " + nBins)
  print("    rmin     =  " + rmin)
  print("    rmrat    =  " + rmrat)
  print("    rho      =  " + rho)
  
  setupbins(nBins, rmin, rmrat, rho, rmass, r, dr)
  
  ; Open the control history file
  history = addfile("/home/bardeen/waccm/carma/analysis/125_lt_gw2_28/h0.2028.mon.nc", "r")
  
  ; Compute the pressure levels ...
  time     = history->time
  nTimes   = dimsizes(time)
  lev     = history->lev
  nLevs   = dimsizes(lev)
  lat     = history->lat
  nLats   = dimsizes(lat)
  lon     = history->lon
  nLons   = dimsizes(lon)
  
  ; Determine an approximate dz/dp.
  hyam 	= history->hyam
  hybm 	= history->hybm
  hyai 	= history->hyai
  hybi 	= history->hybi
  P0 	  = history->P0
  PS 	  = history->PS
  Q     = history->Q
  T     = history->T
  u     = history->U
  v     = history->V
  
  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)"
  
  month = 7 

  ; Calculate the air pressures at the interfaces.
  print("")
  print("  Calculating air pressure and density ...")
  P = Q
  P = (/ pres_hybrid_ccm(PS(:,:,:), P0, hyam(:) ,hybm(:)) /)
  P@units         = "Pa"
  P@long_name     = "pressure"
  
  ; Calculate the density assuming an ideal gas.
  arho            = P
  arho@units      = "g/cm3"
  arho@long_name  = "air density"

  setupair(history, P, arho)
  arho = arho(:,:,:,:) * 1e3 / 1e6
  
  ; Get the Number density (#/cm-3) for each bin at the 
  ; specified time and location.
  print("")
  print("  Calculating dust number, area and mass densities ...")
  N = new((/ nLevs, nBins, nTimes, nLats, nLons /), float)
  N!0 = "lev"
  N&lev = lev
  N!1 = "radius"
  N&radius = r
  N!2 = "time"
  N&time = time
  N!3 = "lat"
  N&lat = lat
  N!4 = "lon"
  N&lon = lon
  N@units         = "#/cm3"
  N@long_name     = "number density"
  getN(history, "DUST", rmass, arho, N)


  ; Calculate the zonal average of total number, area and mass in the
  ; desired bin range.
  firstBin  =  7           ; smallest particle to include in the totals
  lastBin   = nBins-1     ; largest particle to include in the totals
  print("  smallest radius = " + r(firstBin))

  N1 = dim_sum_Wrap(N(time|:, lev|:, lat|:, lon|:, radius|firstBin:lastBin))
  N1@long_name = "Concentration, r~B~0~N~ >= 1 nm"
  N1@units     = "# cm~S~-3~N~"
  
  ; Make polar plots of the evolution of the dust through the polar vortex
  ; in the summer vortex.
  Zs = (/ 91, 88, 85 /)
  nZs = dimsizes(Zs)
  Ps = (/ 1.5381e-3, 2.6173e-3, 4.4568e-3 /)

  hemisphere = "NH"
  if (hemisphere .eq. "NH") then
    minLat     = 45.0
    maxLat     = 90.0
  else
    minLat     = -45.0
    maxLat     = -90.0
  end if

  ; Make a 3 plot panel ...
  print("")
  print("  Plotting ...")
 
  nPlots = nZs
  plot = new((/ nPlots /), graphic)
  
  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnPanelLabelBar     = True
  resP@lbTitleOn            =  True                ; turn on title
  resP@lbTitlePosition      = "Right"              ; title position
  resP@lbTitleFontHeightF   = .015                ; make title smaller
  resP@lbTitleDirection     = "Across"             ; title direction
  resP@lbTitleString        = "cm~S~-3~N~"
;  resP@lbLabelStride        = 2
        
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@gsnPolar              = hemisphere
  if (hemisphere .eq. "NH") then
    res@mpMinLatF             = minLat
  else
    res@mpMaxLatF             = minLat
  end if
  res@mpFillOn              = False
  res@cnFillOn              = True                  ; color fill
  res@cnLinesOn             = False                 ; no contour lines
  res@gsnScalarContour      = True              ; vectors over scalars
  res@vcRefMagnitudeF       = 75.0              ; vec magnitude
  res@vcRefLengthF          = 0.04              ; size of ref vector
  res@vcLineArrowColor      = "grey"           ; white vectors
;  res@vcMinDistanceF        = 0.02              ; thins arrows near pole  
  res@vcMinDistanceF        = 0.035              ; thins arrows near pole  
  res@vcRefAnnoFontHeightF    = 0.025              ; thins arrows near pole  
  res@vcGlyphStyle          = "CurlyVector"     ; turns on curley vectors
;  res@vcRefAnnoArrowLineColor   = "black"       ; black ref vector
;  res@vcRefAnnoArrowUseVecColor = False
  res@vcRefAnnoString1      = "75 (m/s)"
  res@gsnSpreadColors     = True                ; use full colormap
  res@gsnSpreadColorStart = 2 
  res@gsnSpreadColorEnd   = -3 
  res@lbLabelBarOn          = False
  res@cnInfoLabelOn         = False
;  res@gsnLeftStringFontHeightF = 0.025
  res@gsnCenterStringFontHeightF = 0.04
;  res@gsnRightStringFontHeightF = 0.025
  res@gsnPolarLabelFontHeightF = 0.025
  
  months = (/ "Dec", "January", "Feb", "Mar", "April", "May", "June", "July", "Aug", "September", "Oct", "Nov" /)

  ; Open a workstation
  wks = gsn_open_wks("pdf","images/125_polar")
  gsn_define_colormap(wks,"WhBlGrYeRe")                ; choose colormap
  
  print("")
  print("    ND > 1 nm")
  
  res@gsnLeftString      = ""
  res@gsnRightString     = ""
  
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 5.0                   ; contour interval
  res@cnMinLevelValF       = 0.0             ; min level
  res@cnMaxLevelValF       = 100.0           ; max level
      
  do i = 0, nZs-1
    print("  Plotting for " + Zs(i) + " km, " + Ps(i) + " mb ...")
  
    res@gsnCenterString      = Zs(i) + " km"
    plot(i) = gsn_csm_vector_scalar_map_polar(wks, u(month, {Ps(i)}, {minLat:maxLat}, :), v(month, {Ps(i)}, {minLat:maxLat}, :), N1(month, {Ps(i)}, {minLat:maxLat}, :), res)
  end do

  gsn_panel(wks, plot, (/ 1, 3/), resP)
end

monave_za()
;monave_polar()
