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


undef("monave_za")
procedure monave_za(month:integer)

local files, nFiles, historyFile

begin
  print("")
  print("  Plot monthly average dust densities ...")
  
  ; Set up the particle bins ...
  nBins = 21
  rmin  = 1e-7 ; cm
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
  files   = (/ "analysis/mdust/h0.2028.mon.nc", "analysis/mdust_nosed/h0.2028.mon.nc" /)
  nFiles  = dimsizes(files)
  history = addfile(files(0), "r")
  
  lev     = history->lev
  nLevs   = dimsizes(lev)
  lat     = history->lat
  nLats   = dimsizes(lat)
  
  delete(history)
  
  lev_p   = lev
  lev_p@units = "mb"
  lev_p@long_name = "Pressure (mb)"
  
  ; Make three sets of data:
  ;  control
  ;  test
  ;  difference
  runs = (/ "Control", "No Sedimentation", "(No Sedimentation - Control) / Control" /)
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
  N1s@long_name = "Concentration"
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

  
  ; Read in the global average data for one year.
  do i = 0, nFiles - 1
  
    history = addfile(files(i), "r")

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
    firstBin  =  0           ; smallest particle to include in the totals
    lastBin   = nBins-1     ; largest particle to include in the totals
  
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

    Ts(i, :, :) = dim_avg(T(month, :, :, :))
    N1s(i, :, :) = dim_avg(N1(month, :, :, :))
    NDs(i, :, :) = dim_avg(ND(month, :, :, :))
    ADs(i, :, :) = dim_avg(AD(month, :, :, :))
    MDs(i, :, :) = dim_avg(MD(month, :, :, :)) * 1e16
    REs(i, :, :) = dim_avg(RE(month, :, :, :)) * 1000.0

    delete(N1)
    delete(ND)
    delete(AD)
    delete(MD)
    delete(RE)
    delete(N)
  end do
  
  ; Calculate the differences ...
  N1s(nFiles, :, :) = N1s(1, :, :) - N1s(0, :, :)
  NDs(nFiles, :, :) = NDs(1, :, :) - NDs(0, :, :)
  ADs(nFiles, :, :) = ADs(1, :, :) - ADs(0, :, :)
  MDs(nFiles, :, :) = MDs(1, :, :) - MDs(0, :, :)
  REs(nFiles, :, :) = REs(1, :, :) - REs(0, :, :)
  
  ; Check for 90% probability (zc >= 1.645)
  zc = 1.645
  
  ; Read in the 7 year average and standard deviation.
  favg = addfile("analysis/mdust/h0.mon.ave.nc", "r")
  fstd = addfile("analysis/mdust/h0.mon.std.nc", "r")
  
  NDfavg = favg->DUSTND
  NDfstd = fstd->DUSTND

  MDfavg = favg->DUSTMD
  MDfstd = fstd->DUSTMD
  
  ; Calculate a zonal average for the same month.
  ;
  ; NOTE: It is probably better to compute zc on all the data and then
  ; average, but for simpilicity, do a first cut on the averaged data.
  ;
  ; NOTE: Better still is to compute the zonally averaged numbers for
  ; each of the 7 years and then find the standard deviation of that.
  ; However, that would require a fair amount of reprocessing.
  ;
  ; NOTE: Doing it this way probably yields a smaller area > 90% than
  ; you would actually get. If our plots look like they don't have much
  ; significance, then The data should be reprocessed.
  NDavg = dim_avg_Wrap(NDfavg(month, :, :, :))
  NDstd = dim_avg_Wrap(NDfstd(month, :, :, :))
  MDavg = dim_avg_Wrap(MDfavg(month, :, :, :)) * 1e16
  MDstd = dim_avg_Wrap(MDfstd(month, :, :, :)) * 1e16
  
  ; NOTE: If std is 0, then the division will fail. Therefore set
  ; std to a fill value for the calculation. These areas would not
  ; be stippled anyway.
  NDstd@_FillValue = -999.0
  MDstd@_FillValue = -999.0
  NDstd = where(NDstd.eq.0, NDstd@_FillValue, NDstd)
  MDstd = where(MDstd.eq.0, MDstd@_FillValue, MDstd)
  
  ; This sets the dimensions, control variables, ...
  NDzc = NDavg
  NDzc = (/ NDzc * 0.0 /)
  delete(NDzc@long_name)
  delete(NDzc@units)
  MDzc = NDzc
  
  ; Calculate a test zc
  ;
  ; NOTE: std can be 0, so need to break this out ...
  NDzc = abs(NDs(1, :, :) - NDavg(:, :)) / NDstd(:, :) 
  MDzc = abs(MDs(1, :, :) - MDavg(:, :)) / MDstd(:, :)
  

  ; Make a 4 plot panel ...
  print("")
  print("  Plotting ...")
 
  minP      = 5e-5        ; largest pressure to plot ...
  maxP      = 15.0        ; largest pressure to plot ...
  
  nPlots = 4
  plot = new((/ nPlots /), graphic)
  
  resP                      = True
  resP@gsnMaximize          = True
  resP@gsnPanelXWhiteSpacePercent = 5
  
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@lbOrientation         = "Vertical"
  res@pmLabelBarOrthogonalPosF = .1
  res@cnLinesOn             = False
  res@cnFillOn              = True                  ; turn on color fill
  res@cnInfoLabelOn         = False
  res@cnLineLabelsOn        = False
  res@gsnSpreadColors       = True                  ; use full range of colors
  res@cnFixFillBleed        = True
;  res@cnFillMode            = "CellFill"
;  res@lbTitlePosition       = "Right"              ; title position
;  res@lbTitleDirection      = "Across"             ; title direction
  res@lbLabelFontHeightF    = .025                ; make title smaller
  res@lbLeftMarginF         = .15
  
  res@gsnLeftStringFontHeightF = 0.035
;  res@gsnCenterStringFontHeightF = 0.0375
  res@gsnRightStringFontHeightF = 0.035
  res@tiMainFontHeightF     = 0.04
  res@tmXBLabelFontHeightF  = 0.03
  res@tmYLLabelFontHeightF  = 0.03
  res@tmYRLabelFontHeightF  = 0.03
  res@tiYAxisFontHeightF    = 0.035

  res2 = True                            ; res2 probability plots
  res2@gsnDraw             = False       ; Do not draw plot
  res2@gsnFrame            = False       ; Do not advance frome
  res2@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res2@cnMinLevelValF      = 0.00        ; set min contour level
  res2@cnMaxLevelValF      = 2.00        ; set max contour level
  res2@cnLevelSpacingF     = 0.10        ; set contour spacing
  res2@cnInfoLabelOn       = False       ; turn off info label
  res2@cnLinesOn           = False       ; do not draw contour lines
  res2@cnLineLabelsOn      = False       ; do not draw contour labels
  res2@cnFillScaleF        = 0.6         ; add extra density

  months = (/ "Dec", "Jan", "Feb", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov" /)

  ; Open a workstation
  wks = gsn_open_wks("eps","images/paper/2007jd009515r-p10_orig")
  gsn_merge_colormaps(wks, "gui_default", "BlWhRe")
  
  print("")
  print("    number density")
  zonal                     = NDs
  
;  res@lbTitleString        = zonal@units
  res@lbLabelStride         = 1
  res@cnLevelSelectionMode  = "ManualLevels" 

  print("      min :  " + min(zonal(1, {minP:maxP}, :)))
  print("      max :  " + max(zonal(1, {minP:maxP}, :)))

  res@gsnSpreadColorStart = 2 
  res@gsnSpreadColorEnd   = 23 

  res@cnLevelSpacingF       = 500.0                  ; contour interval
  res@cnMinLevelValF        = 500.0                  ; contour interval
  res@cnMaxLevelValF        = 7000.0                 ; contour interval  
  res@tiMainString          = runs(1)
  plot(0) = gsn_csm_pres_hgt(wks, zonal(1, {minP:maxP}, :), res)
 
  print("      min :  " + min(zonal(2, {minP:maxP}, :)))
  print("      max :  " + max(zonal(2, {minP:maxP}, :)))

  res@gsnSpreadColorStart = 24
  res@gsnSpreadColorEnd   = 124 


  res@tiMainString          = runs(2)

  zonal(2, {minP:maxP}, :) = (NDs(2, {minP:maxP}, :) / NDs(0, {minP:maxP}, :)) * 100.

  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/ -55., -45., -35., -25., -15., -5., 5., 15., 25., 35., 45., 55. /)
;  symMinMaxPlt((/ -100, 100/), 11, False, res)
  zonal@units = "%"

  plot(1) = gsn_csm_pres_hgt(wks, zonal(2, {minP:maxP}, :), res)
  plot2   = gsn_csm_contour(wks, NDzc({minP:maxP}, :), res2) 
  plot2   = ShadeLtContour(plot2, zc, 17)  ;
  overlay (plot(1), plot2)
  delete(plot2)
  

  print("")
  print("    mass density")
  zonal = MDs
  zonal&lev = lev_p

;  res@lbTitleString        = zonal@units
  print("      min :  " + min(zonal(1, {minP:maxP}, :)))
  print("      max :  " + max(zonal(1, {minP:maxP}, :)))

  res@gsnSpreadColorStart = 2 
  res@gsnSpreadColorEnd   = 23 

  delete(res@cnLevels)
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/ .2, .4, .8, 1, 2, 4, 8, 10, 20, 40, 80, 100, 200 /)
  res@lbLabelStride        = 1

  res@tiMainString          = ""
  plot(2) = gsn_csm_pres_hgt(wks, zonal(1, {minP:maxP}, :), res)

  print("      min :  " + min(zonal(2, {minP:maxP}, :)))
  print("      max :  " + max(zonal(2, {minP:maxP}, :)))
  delete(res@cnLevels)

  res@gsnSpreadColorStart = 24
  res@gsnSpreadColorEnd   = 124 

;  res@cnLevels = (/ -40, -10, -4, -1, -.4, -.1, 0, .1, .4, 1, 4, 10, 40 /)
;  res@lbLabelStride        = 1
  res@tiMainString          = ""
  
  res@cnLevels = (/ -110, -90, -70, -50, -30, -10, 10, 30, 50, 70, 90, 110 /)
  zonal(2, {minP:maxP}, :) = (MDs(2, {minP:maxP}, :) / MDs(0, {minP:maxP}, :)) * 100.
;  symMinMaxPlt((/-100., 140. /), 11, False, res)
  zonal@units = "%"
 
  plot(3) = gsn_csm_pres_hgt(wks, zonal(2, {minP:maxP}, :), res)
  plot2   = gsn_csm_contour(wks, MDzc({minP:maxP}, :), res2) 
  plot2   = ShadeLtContour(plot2, zc, 17)  ;
  overlay (plot(3), plot2)
  delete(plot2)

  gsn_panel(wks, plot, (/ 2, 2 /), resP)

end

; July
month = 7

monave_za(month)
