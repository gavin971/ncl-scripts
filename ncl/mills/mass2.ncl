; Meteoritic Dust

; Constants
;
; NOTE: These are from CAM, shr_const_mod.F90.
PI      = 4.0 * atan(1.0)         ; pi
G       = 9.80616                 ; acceleration of gravity ~ m/s^2
BOLTZ   = 1.38065e-23             ; Boltzmann's constant ~ J/K/molecule
AVOGAD  = 6.02214e26              ; Avogadro's number ~ molecules/kmole
RGAS    = AVOGAD*BOLTZ            ; Universal gas constant ~ J/K/kmole
MWDAIR  = 28.966                  ; molecular weight dry air ~ kg/kmole
RDAIR   = RGAS/MWDAIR             ; Dry air gas constant ~ J/K/kg
MWH2O   = 18.000                  ; molecular weight of water vapor
PTOP    = 100.0                   ; pressure at the top of the model (non-zero) ~ Pa
REARTH  = 6.37122e6               ; Earth's radius (m)
DG2RAD  = PI / 180.0              ; rad deg^-1

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


undef("calculateMass")
procedure calculateMass(WACCM_RUN:string)

begin
  histPath = "run/" + WACCM_RUN + "/*cam2.h0*.nc"
  files = systemfunc("ls " + histPath)
  nFiles = dimsizes(files)
  
  if (nFiles .le. 0) then
    print("ERROR: No files found for ... " + histPath)
  else 
  
    ; Get the total number of times.
    histories = addfiles(files(:nFiles-4), "r")
    ListSetType(histories, "cat")
    
    time = histories[:]->time
    nTimes = dimsizes(time)
    ltime = 0
    delete(histories)

    ; Get the basic information from the first file.
    history = addfile(files(0), "r")
    
    lev = history->lev
    nLevs = dimsizes(lev)
    lat = history->lat
    nLats = dimsizes(lat)
    lon = history->lon
    nLons = dimsizes(lon)
    
    ; Setup the same bin structure used during the model run.
;    nBins = 21
;    rmin  = 1e-7 ; cm
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

    ; Setup weights for calculating the total mass.
    gw = history->gw

    ; Global ...
    dlon = fabs(lon(1) - lon(0))
    dlat = fabs(lat(1) - lat(0))

    ; We need to do some special handling to calculate these values to
    ; compensate for the poles and the equator.
    dlat2 = dlat / 2.
    dh = new((/ nLats /), typeof(lat))

    do i = 0, nLats - 1

      ; The normal handling, +- dlat / 2 in width.
      sdlat = sin((lat(i) + dlat2) * DG2RAD) - sin((lat(i) - dlat2) * DG2RAD)

      ; At the poles, we only want to dlat lat in one direction.
      if (lat(i) .eq. -90.)  then
        sdlat = 1 + sin((lat(i) + dlat2) * DG2RAD)
      end if

      if (lat(i) .eq. 90.) then
        sdlat =  1 - sin((lat(i) - dlat2) * DG2RAD)
      end if

      ; At the equator, the peak value is at the equator, so 2 times one side.
      if (lat(i) .eq. 0.)
        sdlat = 2. * (sin((lat(i) + dlat2) * DG2RAD))
      end if

      dh(i) = REARTH * sdlat
    end do

    dx = ((2. * PI) / nLons) * REARTH
    dhdx = dh * dx

    wgt = new((/ nLats, nLons /), double)
    wgt = conform(wgt, dhdx, 0)
    
    totalMass = new((/ nBins, nTimes, nLevs /), float)
    totalMass!0 = "radius"
    totalMass&radius = r
    totalMass!1 = "time"
    totalMass&time = time
    totalMass!2 = "lev"
    totalMass&lev = lev
    totalMass@units         = "kg"
    totalMass@long_name     = "mass"

    ; Open the last history files
    print("")
    print("  Processing files ...")
    
    do f = 0, nFiles-2
      print("    " + (/ files(f) /))

      history = addfile(files(f), "r")

      ; Find the last time.
      ftime = history->time
      nFTimes = dimsizes(ftime)
    
      ; Compute the air mass and density.
      ;
      ; NOTE: Convert arho and amass are in cgs units.
      print("      Calculating air pressure and density ...")
      hyai 	= history->hyai
      hybi 	= history->hybi
      P0 	= history->P0
      PS 	= history->PS

      ; Calculate the air pressures at the interfaces.
      P = pres_hybrid_ccm(PS(:,:,:), P0, hyai(:) ,hybi(:))
      delete(PS)
      pdel = P(:, 1:nLevs, :, :) - P(:, 0:nLevs-1, :, :)
      delete(P)

      ; Get the mass (kg) for each bin.
      print("      Calculating dust mass ...")
      q         = history->Q

      mass = new((/ nFTimes, nLevs, nLats, nLons /), float)
      mass!0 = "time"
      mass&time = ftime
      mass!1 = "lev"
      mass&lev = lev
      mass!2 = "lat"
      mass&lat = lat
      mass!3 = "lon"
      mass&lon = lon

      ; NOTE: To keep the memory footprint smaller, calculate the total
      ; mass for each bin one at a time.
      do b = 0, nBins-1
        binName = "DUST" + sprinti("%2.2i", b+1)
;        print("    " + binName)

        ; Get the mass mixing ratio.
        mmr = history->$binName$

        ; Compute the particle concentration.
        ; NOTE: Ignore the mass at the bottom level, since this is were the mass that
        ; sediments will collect. Perhaops we should change the model to get rid of this
        ; mass.
        mass  = (mmr(:, :, :, :) * (1 + q(:, :, :, :)) * (pdel(:, :, :, :) / G))        
        totalMass(b, ltime:ltime+nFTimes-1, :)  = wgt_areasum2(mass(:, :, :, :), wgt, 0)

        delete(mmr)
      end do
      
      ltime = ltime + nFTimes

      delete(mass)
      delete(pdel)
      delete(q)
      delete(ftime)
    end do          

    system("mkdir -p analysis/" + WACCM_RUN)
    system("rm analysis/" + WACCM_RUN + "/mass.nc") 
    outfile         = addfile("analysis/" + WACCM_RUN + "/mass.nc", "c")
    outfile->DUST   = totalMass
  end if
end

undef("plotMass")
procedure plotMass(WACCM_RUN:string)

begin
  ; Plot the mass
  system("mkdir -p images/" + WACCM_RUN)
  history = addfile("analysis/" + WACCM_RUN + "/mass.nc", "r")
  
  lev    = history->lev
  nLevs  = dimsizes(lev)
  r      = history->radius
  nRs    = dimsizes(r)
  time   = history->time
  nTimes = dimsizes(time)
  dust   = history->DUST
  
  ; Calculates some mass totals.
  kNames = (/   " All Altitudes", " < 30km", " > 30km", " < 20km", " 20-50km", " 50-80km", " > 80km" /)
  kLow   = (/       0,      35,     0,       44,        25,        14,       0 /)
  kHigh  = (/ nLevs-2, nLevs-2,    35,  nLevs-2,        44,        25,      14 /)
  nKs    = dimsizes(kNames)
  
  mass         = new((/ nKs, nRs, nTimes /), float)
  mass!0       = "group"
  mass&group   = kNames
  mass!1       = "radius"
  mass&radius  = r 
  mass!2       = "time"
  mass&time    = time
  
  do i = 0, nKs-1
    mass(i, :, :) = dim_sum_Wrap(dust(:, :, kLow(i) : kHigh(i)))
  end do
    
  mass = mass / 1e6
  tMass = dim_sum_Wrap(mass(group|:, time|:, radius|:))
  
  mass@long_name            = "Mass (Mg)"  
  tMass@long_name           = "Mass (Mg)"  
        
  ; Open a workstation
  wks = gsn_open_wks("pdf","images/" + WACCM_RUN + "/mass")
      
  resP                      = True
  resP@gsnMaximize          = True
  resP@txString             = "Meteoric Dust, Total Mass"

  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@vpWidthF              = 0.8
  res@vpHeightF             = 0.4
      
  nYears = 10
  xbValues = new((/ nYears*4+1 /), "integer")
  xbLabels = new((/ nYears*4+1 /), "string")

  do i = 0, nYears-1
    xbValues(i*4+0) = 365*i+0 
    xbValues(i*4+1) = 365*i+90 
    xbValues(i*4+2) = 365*i+181 
    xbValues(i*4+3) = 365*i+273

    xbLabels(i*4+0) = "J" 
    xbLabels(i*4+1) = "A" 
    xbLabels(i*4+2) = "J~C~" + (i+1) 
    xbLabels(i*4+3) = "O" 
  end do

  xbValues(i*4+0) = 365*i+0 
  xbLabels(i*4+0) = "J"
  res@trXMinF     = 0
  res@trXMaxF     = xbValues(i*4+0)
   

  res@tmXBMode              = "Explicit"
  res@tmXBValues            = xbValues
  res@tmXBLabels            = xbLabels
  
  res@trYMinF               = 0
  
  plot = new((/ 3/), graphic)      

  ; Total Mass
  res@trYMaxF               = max(tMass)
  res@tiMainString          = "Total"
  res@tiMainFontHeightF     = 0.02
  
 resB = res
  resC = res
  resD = res
  res@xyLineColors          = (/ "black", "yellow", "blue", "green", "red", "orange", "purple" /)
  res@xyDashPatterns        = (/ 0, 0, 0, 0, 0, 0, 0 /)
  res@xyLineThicknesses     = (/ 4.0, 2.0, 4.0, 2.0, 2.0, 2.0, 2.0 /)
  res@pmLegendDisplayMode   = "Always"
  res@pmLegendSide          = "Right"
  res@pmLegendWidthF        = 0.15
  res@pmLegendHeightF       = 0.3
  res@lgLabelFontHeightF    = 0.02
  res@xyExplicitLegendLabels = kNames

  plot(0) = gsn_csm_xy(wks, time, tMass, res)

  resB@tiMainString          = kNames(0) + ", Per Bin"
  resB@trYMaxF               = max(mass(0, :, :))
  resB@xyLineColors          = (/ "green", "green", "green", "green", "green", "green", "green", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "red", "red", "red", "red", "red", "red", "red" /)
  resB@xyDashPatterns        = (/ 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6 /)
  resB@xyLineThicknesses     = (/ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0 /)

  plot(1) = gsn_csm_xy(wks, time, mass(0, :, :), resB)
      
  resB@tiMainString          = kNames(1) + ", Per Bin"
  resB@trYMaxF               = max(mass(1, :, :))

  plot(2) = gsn_csm_xy(wks, time, mass(1, :, :), resB)
      
  gsn_panel(wks, plot, (/ 3, 1 /), resP)
  
  resC@xyLineColors          = (/ "blue", "purple", "cyan", "black", "green", "orange", "red" /)
  resC@xyDashPatterns        = (/ 0, 0, 0, 0, 0, 0, 0 /)
  resC@xyLineThicknesses     = (/ 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0 /)
  resC@pmLegendDisplayMode   = "Always"
  resC@pmLegendSide          = "Right"
  resC@pmLegendWidthF        = 0.15
  resC@pmLegendHeightF       = 0.3
  resC@lgLabelFontHeightF    = 0.02
  
  resD@xyLineColors          = (/"black", "green", "orange", "red" /)
  resD@xyDashPatterns        = (/ 0, 0, 0, 0 /)
  resD@xyLineThicknesses     = (/ 2.0, 2.0, 2.0, 2.0 /)
  resD@pmLegendDisplayMode   = "Always"
  resD@pmLegendSide          = "Right"
  resD@pmLegendWidthF        = 0.15
  resD@pmLegendHeightF       = 0.3
  resD@lgLabelFontHeightF    = 0.02

  myRadii = floattointeger(r * 1e8) / 10.0

  resC@tiMainString          = kNames(1) + ", Bins 01-07"
  resC@trYMaxF               = 1.1*max(mass(1, 0:6, :))
  resC@xyExplicitLegendLabels = myRadii(0:6)
  plot(0) = gsn_csm_xy(wks, time, mass(1, 0:6, :), resC)

  resC@tiMainString          = kNames(1) + ", Bins 08-14"
  resC@trYMaxF               = 1.1*max(mass(1, 7:13, :))
  resC@xyExplicitLegendLabels = myRadii(7:13)
  plot(1) = gsn_csm_xy(wks, time, mass(1, 7:13, :), resC)
      
  resC@tiMainString          = kNames(1) + ", Bins 15-21"
  resC@trYMaxF               = 1.1*max(mass(1, 14:20, :))
  resC@xyExplicitLegendLabels = myRadii(14:20)
  plot(2) = gsn_csm_xy(wks, time, mass(1, 14:20, :), resC)
      
  gsn_panel(wks, plot, (/ 3, 1 /), resP)
  
  resC@tiMainString          = kNames(1) + ", Bins 04-10"
  resC@trYMaxF               = 1.1*max(mass(1, 3:9, :))
  resC@xyExplicitLegendLabels = myRadii(3:9)
  plot(0) = gsn_csm_xy(wks, time, mass(1, 3:9, :), resC)

  resC@tiMainString          = kNames(1) + ", Bins 11-17"
  resC@trYMaxF               = 1.1*max(mass(1, 10:16, :))
  resC@xyExplicitLegendLabels = myRadii(10:16)
  plot(1) = gsn_csm_xy(wks, time, mass(1, 10:16, :), resC)
      
  resD@tiMainString          = kNames(1) + ", Bins 14-17"
  resD@trYMaxF               = 1.1*max(mass(1, 13:16, :))
  resD@xyExplicitLegendLabels = myRadii(13:16)
  plot(2) = gsn_csm_xy(wks, time, mass(1, 13:16, :), resD)
      
  gsn_panel(wks, plot, (/ 3, 1 /), resP)
  
  
  ; Plot the last few years of data so we can see it better.
  res@trXMinF     = 6*365

  plot(0) = gsn_csm_xy(wks, time, tMass, res)
  
  gsn_panel(wks, plot(0), (/ 1, 1 /), resP)
end


WACCM_RUN = getenv("WACCM_RUN")

if (ismissing(WACCM_RUN)) then
  WACCM_RUN = "results"
end if
print("Using WACCM_RUN: " + WACCM_RUN)
print("")

calculateMass(WACCM_RUN)
;plotMass(WACCM_RUN)
