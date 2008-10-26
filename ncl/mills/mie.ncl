; Polar Mesospheric Clouds - Mie Calculations

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
RHO_I   = 0.93                    ; Ice Density (g/cm-3)

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
    N(:, :, :, :, b) = (mmr * arho) / rmass(b)
    
    delete(mmr)
  end do
end 


undef("calculateTau")
procedure calculateTau(WACCM_RUN: string)

begin
  print("")
  print("  Calculate the PMC optical depth ...")
  print("    year  : " + year)
  print("    month : " + month)


  ; Setup the same bin structure used during the model run.
  nBins = 28
  rmin  = 2e-8 ; cm
  rmrat = 2.6
  rho   = RHO_I  ; g/cm-3
  
  rmass = new(nBins, float)
  r = new(nBins, float)
  dr = new(nBins, float)
  
  print("")
  print("    Assuming ...")
  print("      nBins    =  " + nBins)
  print("      rmin     =  " + rmin)
  print("      rmrat    =  " + rmrat)
  print("      rho      =  " + rho)
  
  setupbins(nBins, rmin, rmrat, rho, rmass, r, dr)


  ; Read in mie values
  ;   radius (um), Qext, Qabs, Qsca, Qback, g, omega
  print("")
  print("  Reading in data from mie_pmc_273_q.txt ...")
  mieData = asciiread("analysis/mie_pmc_273_q.txt", (/ 28, 6 /), "float")
  
  names = (/ "Q~B~ext~N~", "Q~B~abs~N~", "Q~B~sca~N~", "g", "~F33~w~F~~B~o~N~" /)
  nNames = dimsizes(names)
  
  radii = mieData(:, 0) * 1e3
  radii@long_name  = "Radius(nm)"
  nRadii = dimsizes(radii)
  
  print("")
  print("  Calculate the extinction cross section ...")

  ; Calculate extinction cross section (cm-2) ...
  sigma = PI * ((radii(:) * 1e-7)^2) * mieData(:, 1)
  sigma@long_name = "~F33~s~F~~B~ext~N~ (cm~S~2~N~)"


  print("")
  print("  Checking the history files ...")
  histPath = "/ptmp/bardeen/waccm/pmc/run/" + WACCM_RUN + "/*cam2.h2.*.nc"
  files = systemfunc("ls " + histPath)
  nFiles = dimsizes(files)
  
  if (nFiles .le. 0) then
    print("ERROR: No files found for ... " + histPath)
  else 
  
    ; Get the total number of times.
    histories = addfiles(files, "r")
    ListSetType(histories, "cat")
    
    time = addfiles_GetVar(histories, files, "time")
    nTimes = dimsizes(time)
    ltime = 0
    delete(histories)
    
    print("")
    print("  Processing " + nTimes + " times " + " in " + nFiles + " files ...")
    

    ; Get the basic information from the first file.
    history = addfile(files(0), "r")
    
    lev = history->lev
    ilev = history->ilev
    nLevs = dimsizes(lev)
    lat = history->lat
    nLats = dimsizes(lat)
    lon = history->lon
    nLons = dimsizes(lon)
    gw = history->gw
    
    wgt = new((/ nLats, nLons /), double)
    wgt = conform(wgt, gw / (2*nLons), 0)
    
    crice_tod = new((/ nTimes, nLevs, nLats, nLons /), float)
    crice_tod!0 = "time"
    crice_tod&time = time
    crice_tod!1 = "lev"
    crice_tod&lev = lev
    crice_tod!2 = "lat"
    crice_tod&lat = lat
    crice_tod!3 = "lon"
    crice_tod&lon = lon
    crice_tod@units         = ""
    crice_tod@long_name     = "PMC Optical Depth"
    
    crice_tod(:, :, :, :) = 0.0

    crice_ext = crice_tod
    crice_ext@units         = "km-1"
    crice_ext@long_name     = "Volume Extinction Coefficient"
      
    ; Open the history files
    do f = 0, nFiles-1
      print("")
      print("    " + f + " :  " + (/ files(f) /))

      history = addfile(files(f), "r")

      ; Find the last time.
      ftime = history->time
      nFTimes = dimsizes(ftime)
    
      ; Compute the air mass and density.
      ;
      ; NOTE: Convert arho and amass are in cgs units.
      print("      Calculating air pressure and density ...")
      hyam 	= history->hyam
      hybm 	= history->hybm
      P0 	  = history->P0
      PS   	= history->PS
      q     = history->Q
  
      ; Compute the air mass and density.
      ;
      ; NOTE: Convert arho and amass are in cgs units.
      P               = q
      P@units         = "Pa"
      P@long_name     = "Pressure"

      arho            = q
      arho@units      = "g/cm3"
      arho@long_name  = "air density"

      setupair(history, P, arho)
      arho = arho(:,:,:,:) * 1e3 / 1e6
  
      ; Get the Number density (#/cm-3) for each bin at the 
      ; specified time and location.
      print("      Calculating the PMC number density ...")
      ND = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
      getN(history, "CRICE", rmass, arho, ND)
      
      print("      Calculating the PMC optical depth ...")
      ; Get the Z3 at the edges and then DZ in centimeters.
      ;
      ; NOTE: We only care about the pure pressure porition, so just
      ; do a rough conversion from the hybrid coordinate to pressure levels.
      ;
      ; NOTE: This causes a few of the levels at the bottom to have negative
      ; dZ in some layers, since the real P is less that that shown in
      ; ilev.
      Z3  = history->Z3
      Z3i = vinth2p(Z3, hyam, hybm, ilev, PS, 2, P0, 1, True)

      dZ  = (Z3i(:, 0:nLevs-1, :, :) - Z3i(:, 1:nLevs, :, :)) * 100.0
      
      ; NOTE: Some of the lower levels have odd dZs. Just zero them out,
      ; since we don't care about this region anyway.
      dZ(:, nLevs-6:nLevs-1, :, :) = 0.0 
      
      delete(arho)
      delete(P)
      delete(PS)
      delete(Z3)
      delete(Z3i)
      delete(q)
      
      ; Calculate the total optical depth
      do b = 0, nRadii-1
        crice_ext(ltime:ltime+nFTimes-1, :, :, :) = crice_ext(ltime:ltime+nFTimes-1, :, :, :) + sigma(b) * ND(:, :, :, :, b)
      end do
      crice_tod(ltime:ltime+nFTimes-1, :, :, :) = crice_ext(ltime:ltime+nFTimes-1, :, :, :) * dZ(:, :, :, :)
      
      print("        maximum tau in a grid point is " + max(crice_tod(ltime:ltime+nFTimes-1, :, :, :)))

      ltime = ltime + nFTimes

      delete(ND)
      delete(dZ)
      delete(ftime)
    end do          

    ; Calculate a column total optical depth; however, don't include particles below
    ; the PMC range.
    tmp                   = crice_tod(time|:, lat|:, lon|:, lev|:)
    crice_ctod            = dim_sum_Wrap(tmp(:, :, :, {0.0:0.1}))
    crice_ctod@long_name  = "PMC Column Optical Depth"

    ; Convert from cm-1 to km-1
    crice_ext = (/ crice_ext * 1e5 /)    
    print("")
    print("  Writing Results ...")
    ofilename = "analysis/" + WACCM_RUN + "/crice_tau.nc"

    system("rm " + ofilename) 
    outfile               = addfile(ofilename, "c")
    outfile->time         = time
    outfile->gw           = gw
    outfile->P0           = P0
    outfile->hyam         = hyam
    outfile->hybm         = hybm
    outfile->CRICE_EXT    = crice_ext
    outfile->CRICE_TOD    = crice_tod
    outfile->CRICE_CTOD   = crice_ctod
  end if
end


undef("calculateTauDaily")
procedure calculateTauDaily(WACCM_RUN: string, year:integer, month:integer, day:integer)

begin
  print("")
  print("  Calculate the PMC optical depth ...")
  print("    year  : " + year)
  print("    month : " + month)


  ; Setup the same bin structure used during the model run.
  nBins = 28
  rmin  = 2e-8 ; cm
  rmrat = 2.6
  rho   = RHO_I  ; g/cm-3
  
  rmass = new(nBins, float)
  r = new(nBins, float)
  dr = new(nBins, float)
  
  print("")
  print("    Assuming ...")
  print("      nBins    =  " + nBins)
  print("      rmin     =  " + rmin)
  print("      rmrat    =  " + rmrat)
  print("      rho      =  " + rho)
  
  setupbins(nBins, rmin, rmrat, rho, rmass, r, dr)


  ; Read in mie values
  ;   radius (um), Qext, Qabs, Qsca, Qback, g, omega
  print("")
  print("  Reading in data from mie_pmc_273_q.txt ...")
  mieData = asciiread("analysis/mie_pmc_273_q.txt", (/ 28, 6 /), "float")
  
  names = (/ "Q~B~ext~N~", "Q~B~abs~N~", "Q~B~sca~N~", "g", "~F33~w~F~~B~o~N~" /)
  nNames = dimsizes(names)
  
  radii = mieData(:, 0) * 1e3
  radii@long_name  = "Radius(nm)"
  nRadii = dimsizes(radii)
  
  print("")
  print("  Calculate the extinction cross section ...")

  ; Calculate extinction cross section (cm-2) ...
  sigma = PI * ((radii(:) * 1e-7)^2) * mieData(:, 1)
  sigma@long_name = "~F33~s~F~~B~ext~N~ (cm~S~2~N~)"


  print("")
  print("  Checking the history files ...")
;  histPath = "run/" + WACCM_RUN + "/*cam2.h0*.nc"
  histPath = "/ptmp/bardeen/waccm/pmc/run/" + WACCM_RUN + "/*cam2.h2." + year + "-" + sprinti("%0.2i", month) + "-" + day + "*.nc"
;  histPath = "/ptmp/bardeen/waccm/pmc/run/" + WACCM_RUN + "/*cam2.h3." + year + "-" + sprinti("%0.2i", month) + "-" + day + "*.nc"
  files = systemfunc("ls " + histPath)
  nFiles = dimsizes(files)
  
  if (nFiles .le. 0) then
    print("ERROR: No files found for ... " + histPath)
  else 
  
    ; Get the total number of times.
    histories = addfiles(files, "r")
    ListSetType(histories, "cat")
    
    time = addfiles_GetVar(histories, files, "time")
    nTimes = dimsizes(time)
    ltime = 0
    delete(histories)
    
    print("")
    print("  Processing " + nTimes + " times " + " in " + nFiles + " files ...")
    

    ; Get the basic information from the first file.
    history = addfile(files(0), "r")
    
    lev = history->lev
    ilev = history->ilev
    nLevs = dimsizes(lev)
    lat = history->lat
    nLats = dimsizes(lat)
    lon = history->lon
    nLons = dimsizes(lon)
    gw = history->gw
    
    wgt = new((/ nLats, nLons /), double)
    wgt = conform(wgt, gw / (2*nLons), 0)
    
    crice_tod = new((/ nTimes, nLevs, nLats, nLons /), float)
    crice_tod!0 = "time"
    crice_tod&time = time
    crice_tod!1 = "lev"
    crice_tod&lev = lev
    crice_tod!2 = "lat"
    crice_tod&lat = lat
    crice_tod!3 = "lon"
    crice_tod&lon = lon
    crice_tod@units         = ""
    crice_tod@long_name     = "PMC Optical Depth"
    
    crice_tod(:, :, :, :) = 0.0

    crice_ext = crice_tod
    crice_ext@units         = "km-1"
    crice_ext@long_name     = "Volume Extinction Coefficient"
      
    ; Open the history files
    do f = 0, nFiles-1
      print("")
      print("    " + f + " :  " + (/ files(f) /))

      history = addfile(files(f), "r")

      ; Find the last time.
      ftime = history->time
      nFTimes = dimsizes(ftime)
    
      ; Compute the air mass and density.
      ;
      ; NOTE: Convert arho and amass are in cgs units.
      print("      Calculating air pressure and density ...")
      hyam 	= history->hyam
      hybm 	= history->hybm
      P0 	  = history->P0
      PS   	= history->PS
      q     = history->Q
  
      ; Compute the air mass and density.
      ;
      ; NOTE: Convert arho and amass are in cgs units.
      P               = q
      P@units         = "Pa"
      P@long_name     = "Pressure"

      arho            = q
      arho@units      = "g/cm3"
      arho@long_name  = "air density"

      setupair(history, P, arho)
      arho = arho(:,:,:,:) * 1e3 / 1e6
  
      ; Get the Number density (#/cm-3) for each bin at the 
      ; specified time and location.
      print("      Calculating the PMC number density ...")
      ND = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
      getN(history, "CRICE", rmass, arho, ND)
      
      print("      Calculating the PMC optical depth ...")
      ; Get the Z3 at the edges and then DZ in centimeters.
      ;
      ; NOTE: We only care about the pure pressure porition, so just
      ; do a rough conversion from the hybrid coordinate to pressure levels.
      ;
      ; NOTE: This causes a few of the levels at the bottom to have negative
      ; dZ in some layers, since the real P is less that that shown in
      ; ilev.
      Z3  = history->Z3
      Z3i = vinth2p(Z3, hyam, hybm, ilev, PS, 2, P0, 1, True)

      dZ  = (Z3i(:, 0:nLevs-1, :, :) - Z3i(:, 1:nLevs, :, :)) * 100.0
      
      ; NOTE: Some of the lower levels have odd dZs. Just zero them out,
      ; since we don't care about this region anyway.
      dZ(:, nLevs-6:nLevs-1, :, :) = 0.0 
      
      delete(arho)
      delete(P)
      delete(PS)
      delete(Z3)
      delete(Z3i)
      delete(q)
      
      ; Calculate the total optical depth
      do b = 0, nRadii-1
        crice_ext(ltime:ltime+nFTimes-1, :, :, :) = crice_ext(ltime:ltime+nFTimes-1, :, :, :) + sigma(b) * ND(:, :, :, :, b)
      end do
      crice_tod(ltime:ltime+nFTimes-1, :, :, :) = crice_ext(ltime:ltime+nFTimes-1, :, :, :) * dZ(:, :, :, :)
      
      print("        maximum tau in a grid point is " + max(crice_tod(ltime:ltime+nFTimes-1, :, :, :)))

      ltime = ltime + nFTimes

      delete(ND)
      delete(dZ)
      delete(ftime)
    end do          

    ; Calculate a column total optical depth; however, don't include particles below
    ; the PMC range.
    tmp                   = crice_tod(time|:, lat|:, lon|:, lev|:)
    crice_ctod            = dim_sum_Wrap(tmp(:, :, :, {0.0:0.1}))
    crice_ctod@long_name  = "PMC Column Optical Depth"

    ; Convert from cm-1 to km-1
    crice_ext = (/ crice_ext * 1e5 /)    
    print("")
    print("  Writing Results ...")
    ofilename = "analysis/" + WACCM_RUN + "/crice_tau_" + year + "-" + sprinti("%0.2i", month) + "-" + day + "0.nc"

    system("rm " + ofilename) 
    outfile               = addfile(ofilename, "c")
    outfile->time         = time
    outfile->gw           = gw
    outfile->P0           = P0
    outfile->hyam         = hyam
    outfile->hybm         = hybm
    outfile->CRICE_EXT    = crice_ext
    outfile->CRICE_TOD    = crice_tod
    outfile->CRICE_CTOD   = crice_ctod
  end if
end


undef("plotPhases")
procedure plotPhases()

begin
  print("")
  print("")
  print("  Plot Mie phase functions (normalized) at 273 nm ...")

  ; Read in mie values
  ;   radius (um), Qext, Qabs, Qsca, Qback, g, omega
  print("")
  print("  Read in data from mie_pmc_273_q.txt ...")
  mieData = asciiread("analysis/mie_pmc_273_q.txt", (/ 28, 6 /), "float")
  miePhase = asciiread("analysis/mie_pmc_273_pn.txt", (/ 28, 181 /), "float")
  
  radii = mieData(:, 0) * 1e3
  radii@long_name  = "Radius(nm)"
  nRadii = dimsizes(radii)
  
  miePhase@long_name = "P(~F33~Q~F~)"
  
  angle = ispan(0, 180,1)
  angle@long_name = "~F33~Q~F~ (degrees)"
  
  ; Open a workstation
  wks = gsn_open_wks("pdf","images/pmc2/mie_phase")
  
  hPlot = new((/ nRadii /), graphic)
      
  resP                      = True
  resP@gsnMaximize          = True
  resP@txString             = "Phase Function for PMCs, ~F33~l~F~=273 nm"
  resP@gsnPanelXWhiteSpacePercent = 3
  resP@gsnPanelYWhiteSpacePercent = 3
  
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@trXMinF               = 0
  res@trXMaxF               = 180
  res@xyLineThicknessF     = 4.0
  res@tiMainFontHeightF     = 0.045
  res@tmXBLabelFontHeightF  = 0.04
  res@tmYLLabelFontHeightF  = 0.04
  res@tiYAxisFontHeightF    = 0.05
  res@tiXAxisFontHeightF    = 0.05
  
  do i = 0, nRadii-1    
    res@tiMainString = "" + radii(i) + " nm"
    
    hPlot(i) = gsn_csm_xy(wks, angle, miePhase(i, :), res)
  end do

  gsn_panel(wks, hPlot, (/ 4, 7 /), resP)  
  delete(wks)
end


undef("plotScattering")
procedure plotScattering()

begin
  print("")
  print("")
  print("  Plot Mie phase functions (unnormalized) at 273 nm ...")

  ; Read in mie values
  ;   radius (um), Qext, Qabs, Qsca, Qback, g, omega
  print("")
  print("  Read in data from mie_pmc_273_q.txt ...")
  mieData = asciiread("analysis/mie_pmc_273_q.txt", (/ 28, 6 /), "float")
  miePhase = asciiread("analysis/mie_pmc_273_p.txt", (/ 28, 181 /), "float")
  
  radii = mieData(:, 0) * 1e3
  radii@long_name  = "Radius(nm)"
  nRadii = dimsizes(radii)
  
  miePhase@long_name = "F(~F33~Q~F~)"
  
  angle = ispan(0, 180,1)
  angle@long_name = "~F33~Q~F~ (degrees)"
  
  ; Open a workstation
  wks = gsn_open_wks("pdf","images/pmc2/mie_scattering")
  
  hPlot = new((/ nRadii /), graphic)
      
  resP                      = True
  resP@gsnMaximize          = True
  resP@txString             = "Scattering Function for PMCs, ~F33~l~F~=273 nm"
  resP@gsnPanelXWhiteSpacePercent = 3
  resP@gsnPanelYWhiteSpacePercent = 3
  
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@trXMinF               = 0
  res@trXMaxF               = 180
  res@xyLineThicknessF     = 4.0
  res@tiMainFontHeightF     = 0.045
  res@tmXBLabelFontHeightF  = 0.04
  res@tmYLLabelFontHeightF  = 0.04
  res@tiYAxisFontHeightF    = 0.05
  res@tiXAxisFontHeightF    = 0.05
  
  do i = 0, nRadii-1    
    res@tiMainString = "" + radii(i) + " nm"
    
    hPlot(i) = gsn_csm_xy(wks, angle, miePhase(i, :), res)
  end do

  gsn_panel(wks, hPlot, (/ 4, 7 /), resP)  
  delete(wks)
end


undef("plotProps")
procedure plotProps()

begin
  print("")
  print("")
  print("  Plot Mie properties at 273 nm ...")

  ; Read in mie values
  ;   radius (um), Qext, Qabs, Qsca, Qback, g, omega
  print("")
  print("  Reading in data from mie_pmc_273_q.txt ...")
  mieData = asciiread("analysis/mie_pmc_273_q.txt", (/ 28, 6 /), "float")
  
  names = (/ "Q~B~ext~N~", "Q~B~abs~N~", "Q~B~sca~N~", "g", "~F33~w~F~~B~o~N~" /)
  nNames = dimsizes(names)
  

  ; Open a workstation
  wks = gsn_open_wks("pdf","images/pmc2/mie_props")
  
  hPlot = new((/ nNames+1 /), graphic)
      
  resP                      = True
  resP@gsnMaximize          = True
  resP@txString             = "Mie Properties for PMCs, ~F33~l~F~=273 nm"
  resP@gsnPanelXWhiteSpacePercent = 7
  resP@gsnPanelYWhiteSpacePercent = 5
  
  res                       = True
  res@gsnDraw               = False
  res@gsnFrame              = False
;  res@trYMinF               = 0.0
;  res@trYMaxF               = 4.0
  res@trXMinF               = 2e-1
  res@tmXBMinorPerMajor     = 8
  res@tmXTMinorPerMajor     = 8
  res@trXLog                = True
  res@xyLineThicknessF      = 4.0
  
  radii = mieData(:, 0) * 1e3
  radii@long_name  = "Radius(nm)"
  
  do i = 0, nNames-1
    data = mieData(:, i+1)
    data@long_name = names(i)
    
    if (i.eq.1) then
      res@tmYLFormat            = "0@*+^e"
    else
      res@tmYLFormat            = "0@*+^sg"
    end if

    if (i.ge.nNames-2) then
      res@trYMinF               = -1.0
      res@trYMaxF               = 1.0
    end if

    if (i.ge.nNames-1) then
      res@trYMinF               = 0.0
      res@trYMaxF               = 1.0
    end if
    
    hPlot(i) = gsn_csm_xy(wks, radii, data, res)
  end do
  
  ; Calculate extinction cross section ...
  sigma = PI * ((radii * 1e-7)^2) * mieData(:, 1)
  sigma@long_name = "~F33~s~F~~B~ext~N~ (cm~S~2~N~)"
  
  delete(res@trYMinF)
  delete(res@trYMaxF)
  
  res@trYLog                = True
  res@tmYLFormat            = "0@*+^e"
  
  hPlot(i) = gsn_csm_xy(wks, radii, sigma, res)

  gsn_panel(wks, hPlot, (/ 2, 3 /), resP)
  delete(wks)
end



; Main
WACCM_RUN = getenv("WACCM_RUN")

if (ismissing(WACCM_RUN)) then
  WACCM_RUN = "results"
end if
print("Using WACCM_RUN: " + WACCM_RUN)
print("")


calculateTau(WACCM_RUN)

;do month = 5, 9
;  do day = 0, 3
;    calculateTau(WACCM_RUN, 1995, month, day)
;  end do  
;end do

;plotPhases()
;plotScattering()
;plotProps()

