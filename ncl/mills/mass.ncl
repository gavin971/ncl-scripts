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

load "tests/core.ncl"

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
procedure calculateMass()

begin
  histPath = "run/mdust/*cam2.h0*.nc"
  files = systemfunc("ls " + histPath)
  nFiles = dimsizes(files)
  
  if (nFiles .le. 0) then
    print("ERROR: No files found for ... " + histPath)
  else 
    
    ; Open the last history files
    print("  Opening files ...")
    print("    " + (/ files(:nFiles-2) /))

    history = addfiles(files(:nFiles-2), "r")

    ; Find the last time.
    time = history[:]->time
    nTimes = dimsizes(time)
    lev = history[0]->lev
    nLevs = dimsizes(lev)
    lat = history[0]->lat
    nLats = dimsizes(lat)
    lon = history[0]->lon
    nLons = dimsizes(lon)

    ; Setup the same bin structure used during the model run.
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

    ; Setup weights for calculating the total mass.
    ; 
    ; NOTE: For now, just do the last timestep.
    gw = history[0]->gw

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
    
    ; Compute the air mass and density.
    ;
    ; NOTE: Convert arho and amass are in cgs units.
    print("")
    print("  Calculating air pressure and density ...")
    hyai 	= history[0]->hyai
    hybi 	= history[0]->hybi
    P0 	= history[0]->P0
    PS 	= history[:]->PS

    ; Calculate the air pressures at the interfaces.
    P = pres_hybrid_ccm(PS(:,:,:), P0, hyai(:) ,hybi(:))
    delete(PS)
    pdel = P(:, 1:nLevs, :, :) - P(:, 0:nLevs-1, :, :)
    delete(P)

    ; Get the mass (kg) for each bin.
    print("") 
    print("  Calculating dust mass ...")
    q         = history[:]->Q
    mass = new((/ nTimes, nLevs, nLats, nLons /), float)
    mass!0 = "time"
    mass&time = time
    mass!1 = "lev"
    mass&lev = lev
    mass!2 = "lat"
    mass&lat = lat
    mass!3 = "lon"
    mass&lon = lon

    totalMass = new((/ nBins, nTimes /), float)
    totalMass!0 = "radius"
    totalMass&radius = r
    totalMass!1 = "time"
    totalMass&time = time
    totalMass@units         = "kg"
    totalMass@long_name     = "mass"

    ; NOTE: To keep the memory footprint smaller, calculate the total
    ; mass for each bin one at a time.
    do b = 0, nBins-1
      binName = "DUST" + sprinti("%2.2i", b+1)
      print("    " + binName)

      ; Get the mass mixing ratio.
      mmr = history[:]->$binName$

      ; Compute the particle concentration.
      ; NOTE: Ignore the mass at the bottom level, since this is were the mass that
      ; sediments will collect. Perhaops we should change the model to get rid of this
      ; mass.
      mass             = (mmr(:, :, :, :) * (1 + q(:, :, :, :)) * (pdel(:, :, :, :) / G))        
      columnMass       = dim_sum_Wrap(mass(time|:, lat|:, lon|:, lev|0:nLevs-2))       
      totalMass(b, :)  = wgt_areasum2(columnMass(:, :, :), wgt, 0)

      delete(mmr)
      delete(columnMass)
    end do

    delete(mass)
    delete(pdel)
    delete(q)           

    system("rm analysis/mass.nc") 
    outfile         = addfile("analysis/mass.nc", "c")
    outfile->DUST   = totalMass
  end if
end

undef("plotMass")
procedure plotMass()

begin
  ; Plot the mass
  history = addfile("analysis/mass.nc", "r")
  
  r    = history->radius
  time = history->time
  mass = history->DUST
  
  mass = mass / 1e6
  mass@long_name            = "Mass (Mg)"  
        
  ; Open a workstation
  wks = gsn_open_wks("pdf","images/mass")
      
  resP                      = True
  resP@gsnMaximize          = True
  resP@txString             = "Mass of Meteoric Dust"

  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@vpWidthF              = 0.8
  res@vpHeightF             = 0.4
      
  nYears = 10
  xbValues = new((/ nYears*4 /), "integer")
  xbLabels = new((/ nYears*4 /), "string")

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

  res@tmXBMode              = "Explicit"
  res@tmXBValues            = xbValues
  res@tmXBLabels            = xbLabels
  
  res@trYMinF               = 0
  
  plot = new((/ 2 /), graphic)      

  tm = dim_sum_Wrap(mass(time|:, radius|:))
  tm@units                  = "Mg"
  tm@long_name              = "Mass (Mg)"
  res@trYMaxF               = max(tm)
  res@tiMainString          = "Total"
  res@tiMainFontHeightF     = 0.02

  plot(0) = gsn_csm_xy(wks, time, tm, res)

  res@tiMainString          = "Per Bin"
  res@trYMaxF               = max(mass)
  res@xyLineColors          = (/ "green", "green", "green", "green", "green", "green", "green", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "red", "red", "red", "red", "red", "red", "red" /)
  res@xyDashPatterns        = (/ 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6 /)
;      res@xyLineThicknesses     = (/ 4.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0 /)

  plot(1) = gsn_csm_xy(wks, time, mass, res)
      
  gsn_panel(wks, plot, (/ 2, 1 /), resP)
end

calculateMass()
plotMass()
