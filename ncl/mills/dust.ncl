; Meteoritic Dust

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

load "tests/core.ncl"

; Determine the mass in each bin based upon a lognormal size
; distribution.
;
; NOTE: This logic is from CARMA initaer.f.

undef("initaer")
function initaer(rmass[*]:float, r[*]:float, dr[*]:float, n:float, r0:float, rsig:float)

local arg1, arg2, i, nbin, totn

begin

  ; Determine the number of bins.
  nbin = dimsizes(r)
  
  mass = new(nbin, double)
  
;  total = 0.
;  do i = 0, nbin-1
;    arg1 = dr(i) / (sqrt(2.*PI) * r(i) * log(rsig))
;    arg2 = -((log(r(i) / r0))^2) / (2.*(log(rsig))^2)
;    total  = total + arg1 * exp(arg2)
;  end do
;  totn = n / total
totn = n

  do i = 0, nbin-1
    arg1 = totn * dr(i) / (sqrt(2.*PI) * r(i) * log(rsig))
    arg2 = -((log(r(i) / r0))^2) / (2.*(log(rsig))^2)
    
    mass(i) = arg1 * exp(arg2) * rmass(i)
  end do
  
  return(mass)
end


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


; Add some analysis fields to the output files. These include:
;
;
undef("analyzeHistory")
procedure analyzeHistory(histPath:string)

local files, nFiles, historyFile

begin
  files = systemfunc("ls " + histPath)
  nFiles = dimsizes(files)
  
  if (nFiles .le. 0) then
    print("ERROR: No files found for ... " + histPath)
  else 
    
      ; Open the last history files
      print("  Opening file ...")
      print("    " + (/ files(nFiles-1) /))
    
      history = addfile(files(nFiles-1), "r")
      
      ; Find the last time.
      q = history->Q
      Z = history->Z3
      time = history->time
      nTimes = dimsizes(time)
      lev = history->lev
      nLevs = dimsizes(lev)
      lat = history->lat
      nLats = dimsizes(lat)
      lon = history->lon
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
      print("    Assuming ...")
      print("      nBins    =  " + nBins)
      print("      rmin     =  " + rmin)
      print("      rmrat    =  " + rmrat)
      print("      rho      =  " + rho)
      
      setupbins(nBins, rmin, rmrat, rho, rmass, r, dr)
      
      ; Compute the air mass and density.
      ;
      ; NOTE: Convert arho and amass are in cgs units.
      print("    Calculating air pressure and density ...")
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
      print("    Calculating dust number density ...")
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
      
      ; Create averages.
      ; 
      ; NOTE: For now, just do the last timestep.
      print("      Calculating averages ...")
      gw = history->gw
      
      ; Global ...
      print("        global ...")
      Nga = wgt_areaave(N(:, :, nTimes-1, :, :), gw(:), 1.0, 0)        
      Zga = wgt_areaave(Z(nTimes-1, :, :, :), gw(:), 1.0, 0)        
      Pga = wgt_areaave(P(nTimes-1, :, :, :), gw(:), 1.0, 0)        
      
      ; Polar Averages ...
      print("        polar ...")
      Nnp = wgt_areaave(N(:, :, nTimes-1, {60.0:90.0}, :), gw({60.0:90.0}), 1.0, 0)        
      Znp = wgt_areaave(Z(nTimes-1, :, {60.0:90.0}, :), gw({60.0:90.0}), 1.0, 0)        
      Nsp = wgt_areaave(N(:, :, nTimes-1, {-90.0:-60.0}, :), gw({-90.0:-60.0}), 1.0, 0)        
      Zsp = wgt_areaave(Z(nTimes-1, :, {-90.0:-60.0}, :), gw({-90.0:-60.0}), 1.0, 0)
      
      
      ; Scale to dN/dlogr
      do b = 0, nBins-1
        Nga(:, b) = (Nga(:,b) / dr(b)) * r(b) * log10(exp(1.0))
        Nnp(:, b) = (Nnp(:,b) / dr(b)) * r(b) * log10(exp(1.0))
        Nsp(:, b) = (Nsp(:,b) / dr(b)) * r(b) * log10(exp(1.0))
      end do

      
      ; Create variables for the plots.
      ;
      ; NOTE: This is probably not a good choice for the pressure.
      print("      Plotting ...")    
      lev_p = lev
      lev_p@units = "hecto-pascals"
      lev_p@long_name = "Pressure (hPa)"
  
      Nga!0 = "lev_p"
      Nga&lev_p = lev_p
      Nga!1 = "radius"
      Nga&radius = r
      Nga@units = "# cm-3"
      Nga@long_name = "dN/dlogr (# cm-3)"
      
      Nnp!0 = "lev_p"
      Nnp&lev_p = lev_p
      Nnp!1 = "radius"
      Nnp&radius = r
      Nnp@units = "# cm-3"
      Nnp@long_name = "dN/dlogr (# cm-3)"

      Nsp!0 = "lev_p"
      Nsp&lev_p = lev_p
      Nsp!1 = "radius"
      Nsp&radius = r
      Nsp@units = "# cm-3"
      Nsp@long_name = "dN/dlogr (# cm-3)"
      
      ; Open a workstation
      wks = gsn_open_wks("pdf","images/dxdlogr")
      
      resP                      = True
      resP@gsnMaximize          = True
      resP@gsnPaperOrientation  = "portrait"
      resP@txString             = "Meteoric Dust, Global, NP & SP Averages"

      res                       = True
      res@gsnDraw               = False
      res@gsnFrame              = False
      res@txFuncCode            = "~"
      res@trYLog                = True
      res@trXLog                = True
      res@tmXBMinorPerMajor     = 8
      res@trXMinF               = 0.9
      res@trXMaxF               = 110.0
      
      rum = r * 1e7
      rum@units = "nm"
      rum@long_name = "Radius (nm)"

      plot = new((/ 9 /), graphic)
      dataN = new((/ 3, nLevs, nBins /), float)
      dataN(0, :, :) = Nga
      dataN(1, :, :) = Nnp
      dataN(2, :, :) = Nsp
      
      ; For log plots, we need the minimum value that is greater than 0.
;      d1d = ndtooned(dataN)
;      res@trYMinF               = min(d1d(ind(d1d .gt. 0.0)))
;      delete(d1d)
      res@trYMaxF               = max(dataN)
      res@trYMinF               = res@trYMaxF / 1e10
      res@xyLineColors          = (/ "black", "blue", "red" /)
      res@xyDashPatterns        = (/ 0, 0, 0 /)
      res@xyLineThicknesses     = (/ 2.0, 2.0, 2.0 /)
      
      ks = (/ 11, 14, 18, 21, 25, 30, 35, 44, 54/)
      nKs = dimsizes(ks)
      
      do i = 0, nKs-1
        k = ks(i)
        res@tiMainString = "P=" + Pga(k)*1e-2 + "(mb), Z=" + Zga(k)*1e-3 + "(km)"
        plot(i) = gsn_csm_xy(wks, rum, dataN(:, k, :), res)
      end do
      
      gsn_panel(wks, (/ plot(0), plot(1), plot(2), plot(3), plot(4), plot(5), plot(6), plot(7), plot(8) /), (/ 3, 3 /), resP)

      ; Now plot surface area density (um2 cm-3)
      dataA = dataN
      do b = 0, nBins-1
       dataA(:, :, b) = dataA(:, :, b) * 4. * PI * (r(b)^2) * 1e8
      end do
      
;      d1d = ndtooned(dataA)
;      res@trYMinF               = min(d1d(ind(d1d .gt. 0.0)))
;      delete(d1d)
      res@trYMaxF               = max(dataA)
      res@trYMinF               = res@trYMaxF / 1e10
      dataA@units                = "um2 cm-3"
      dataA@long_name            = "dA/dlogr (um2 cm-3)" 

      do i = 0, nKs-1
        k = ks(i)
        res@tiMainString = "P=" + Pga(k)*1e-2 + "(mb), Z=" + Zga(k)*1e-3 + "(km)"
        plot(i) = gsn_csm_xy(wks, rum, dataA(:, k, :), res)
      end do

      gsn_panel(wks, (/ plot(0), plot(1), plot(2), plot(3), plot(4), plot(5), plot(6), plot(7), plot(8) /), (/ 3, 3 /), resP)

      ; Now plot mass density (um3 cm-3)
      dataM = dataN
      do b = 0, nBins-1
       dataM(:, :, b) = dataM(:, :, b) * rmass(b)
      end do
      
;      d1d = ndtooned(dataM)
;      res@trYMinF               = min(d1d(ind(d1d .gt. 0.0)))
;      delete(d1d)
      res@trYMaxF               = max(dataM)
      res@trYMinF               = res@trYMaxF / 1e10
      dataM@units                = "g cm-3"
      dataM@long_name            = "dM/dlogr (g cm-3)" 

      do i = 0, nKs-1
        k = ks(i)
        res@tiMainString = "P=" + Pga(k)*1e-2 + "(mb), Z=" + Zga(k)*1e-3 + "(km)"
        plot(i) = gsn_csm_xy(wks, rum, dataM(:, k, :), res)
      end do

      gsn_panel(wks, (/ plot(0), plot(1), plot(2), plot(3), plot(4), plot(5), plot(6), plot(7), plot(8) /), (/ 3, 3 /), resP)
  end if
end


; Add tracers to the initial conditions file.
;
; NOTE: The mass is defined as in g/cm3.
;
; NOTE: Fill the troposphere & stratosphere ~1-1000 mb.

undef("addtracer_particle")
procedure addtracer_particle(icPath:string, shortName:string, mixtype:string, mass[*]:double)

local i, j, k, b, pdel, pdeldry

begin

  ; Determine the number of bins.
  nbin = dimsizes(mass)
  
  ; Open the initial conditions file.
  icFile = addfile(icPath, "w")
  
  ; Read in the dimensions.
  lev 	= icFile->lev
  lat 	= icFile->lat
  lon 	= icFile->lon
  
  nlev = dimsizes(lev)
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  
  ; Do we need to iterate (or do array math) to get ptdeldry at each
  ; grid location and vertical level.
  hyai 	= icFile->hyai
  hybi 	= icFile->hybi
  P0 	= icFile->P0
  PS 	= icFile->PS
  Q     = icFile->Q
  T     = icFile->T
 
  ; Calculate the air pressures at the interfaces.
  pint = pres_hybrid_ccm(PS(:,:,:), P0, hyai(:) ,hybi(:))

  ; Determine the pressures and vertical scale. The vertical spacing
  ; (dz) is based upon hydrostatic balance.
  ;
  ; NOTE: Some models use 0 as the top level, so limit the top pressure
  ; to PTOP in the dz calculation
  pdel    = new((/ nlev, nlat, nlon /), double)
  pdeldry = new((/ nlev, nlat, nlon /), double)
  dz      = new((/ nlev, nlat, nlon /), double)
  
  do i = 0, nlon-1
    do j = 0, nlat-1
      do k = 0, nlev-1        
        pdel(k, j, i)    = pint(0, k+1, j, i) - pint(0, k, j, i)
        dz(k, j, i) = - (RDAIR * T(0, k, j, i) / G) * log(max((/ PTOP, pint(0, k, j, i) /) / pint(0, k+1, j, i)))
      end do
    end do
  end do
  
  pdeldry(:, :, :) = pdel(:, :, :) * (1. - Q(time|0, lev|:, lat|:, lon|:))

  ; Determine the mass of air (kg m-2).
  if (mixtype .eq. "dry") then
    amass = pdeldry / G
  else
    amass = pdel / G
  end if

  ; Need to iterate over each bin.
  tmp = new((/ 1, nlev, nlat, nlon /), double)
  tmp!0 = "time"
  tmp!1 = "lev"
  tmp&lev = lev
  tmp!2 = "lat"
  tmp&lat = lat
  tmp!3 = "lon"
  tmp&lon = lon

  dz!0 = "lev"
  dz&lev = lev
  dz!1 = "lat"
  dz&lat = lat
  dz!2 = "lon"
  dz&lon = lon
  
  amass!0 = "lev"
  amass&lev = lev
  amass!1 = "lat"
  amass&lat = lat
  amass!2 = "lon"
  amass&lon = lon
  
  ; Alternate 
  ptop = (/ 100.,   1.,  .01,   1.,    .01,    .01 /)
  pbot = (/ 992., 100., 1.  , 992., 100.  , 992.   /)
  
  do b = 0, nbin-1

   n = b % 6
    
    ; Name the bin.
    name = shortName + sprinti("%2.2i", b+1)
  
    ; Create the field.
    filevardef(icFile, name, "double", (/"time","lev","lat","lon"/))
    icFile->$name$@units        = "kg/kg"
    icFile->$name$@long_name    = "CARMA particle " + name

    ; Determine the mass mixing ratio.
    ;
    ; NOTE: Convert mass from g cm-3 to kg m-2, since amass
    ; is in kg m-2.
    ;
    ; mmr (kg/kg) = mass (g/cm3) * 1e-3 (kg/g) * 1e4 (cm2/m2) * dz (m) * 100 (cm/m) / amass (kg/m2)
    tmp = 0.0
    tmp(0, {ptop(n):pbot(n)}, :, :) = (mass(b) * 1.e3 * dz({ptop(n):pbot(n)}, :, :)) / amass({ptop(n):pbot(n)}, :, :)

    ; Write the field.
    icFile->$name$ = (/ tmp /)
  end do
  
  delete(tmp)
end


; Generate the initial tracer distribution.
;
; NOTE: The values for nbin, rmin, rmrat. and rho should match those
; in CARMA (aerrad.h and defineaer.f).
undef("newic")
procedure newic(srcIcPath:string)

begin
  system("mkdir -p dyn")

  icPath    = "dyn/ic.nc"
  
  ; Copy the initial conditions file.
  print("  Creating " + icPath)
  system("cp " + srcIcPath + " " + icPath)
  
  ; Add initial conditions for Sulfate CN from a
  ; lognormal distribution.
  ;
  ; NOTE: Only add these values in the tropopause region
  ; This will eliminate cloud formation in other areas.
  shortName = "DUST"
  mixtype   = "wet"
  nbin      = 21
  rmin      = 1e-7              ; cm
  rmrat     = 2.0
  rho       = 2.0               ; g/cm3
    
  print("  Adding particles ... " + shortName)
  rmass = new(nbin, float)
  r = new(nbin, float)
  dr = new(nbin, float)
  setupbins(nbin, rmin, rmrat, rho, rmass, r, dr)

  ; Put the same number of particles in each bin
;  N = 10                       ; #/cm3
;  mass = new(nbin, double)
  
;  mass(:) = N * rmass(:)
;  addtracer_particle(icPath, shortName, mixtype, mass)
  
  ; Put the same mmr in, which simulates a well mixed species.
  mmr = 1e-13
    
  ; Open the initial conditions file.
  icFile = addfile(icPath, "w")

  ; Read in the dimensions.
  lev 	= icFile->lev
  lat 	= icFile->lat
  lon 	= icFile->lon
  
  nlev = dimsizes(lev)
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)

  tmp = new((/ 1, nlev, nlat, nlon /), double)
  tmp!0 = "time"
  tmp!1 = "lev"
  tmp&lev = lev
  tmp!2 = "lat"
  tmp&lat = lat
  tmp!3 = "lon"
  tmp&lon = lon

  ; However, the alternate the pattern with each bin.
  ptop = (/ 100.,   1., .01, 6e-6,   1.,  .01, 6e-6,  .01, 6e-6, 6e-6 /)
  pbot = (/ 992., 100.,  1.,  .01, 992., 100.,   1., 992., 100., 992. /)
  
  do b = 0, nbin-1

    n = b % 10
    
    ; Name the bin.
    name = shortName + sprinti("%2.2i", b+1)
  
    ; Create the field.
    filevardef(icFile, name, "double", (/"time","lev","lat","lon"/))
    icFile->$name$@units        = "kg/kg"
    icFile->$name$@long_name    = "CARMA particle " + name

    ; Set the mmr.
    tmp = 0.0
    tmp(0, {ptop(n):pbot(n)}, :, :) = mmr

    ; Write the field.
    icFile->$name$ = (/ tmp /)
  end do
end

