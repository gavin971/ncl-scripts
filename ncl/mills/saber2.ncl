; Saber Temperatures

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
    N(:, :, :, :, b) = (mmr * arho) / rmass(b)
    
    delete(mmr)
  end do
end 


undef("gridSaber")
procedure gridSaber(gridFile:string)

begin
  print("")
  print("  Grid Saber data daily ...")

  print("")
  print("  Setup Grid ...")
  ; Get the latitude and longitude grids.
  history = addfile(gridFile, "r")
  
  startDay = day_of_year(2007, 5, 1)
  stopDay = day_of_year(2007, 7, 31)
  nTimes =  stopDay - startDay + 1
  time = fspan(startDay, stopDay, nTimes)
  time@long_name = "time"
  time@units = "days since 2007-01-01 00:00:00"
  
  nLevs = 151
  lev = fspan(0., 150., nLevs)

  lat = history->lat
  nLats = dimsizes(lat)
  lon = history->lon
  nLons = dimsizes(lon)
  dlon = abs(lon(0) - lon(1))
  
  gw = history->gw
  
  ; Setup the output variables
  nPoints = new((/ nTimes, nLevs, nLats, nLons /), float)
  nPoints!0 = "time"
  nPoints&time = time
  nPoints!1 = "z"
  nPoints&z = lev
  nPoints!2 = "lat"
  nPoints&lat = lat
  nPoints!3 = "lon"
  nPoints&lon = lon
  nPoints@long_name = "Number of Observations"
  nPoints@_FillValue = -999

  nPoints(:, :, :, :) = 0
  
  T = new((/ nTimes, nLevs, nLats, nLons /), float)
  T!0 = "time"
  T&time = time
  T!1 = "z"
  T&z = lev
  T!2 = "lat"
  T&lat = lat
  T!3 = "lon"
  T&lon = lon
  T@long_name = "Temperature"
  T@units = "K"
  T@_FillValue = -999

  T(:, :, :, :) = 0.
  
    ; Process the data one day at a time
  print("")
  print("  Processing files ...")

  histPath = "~/Data/SABER/*20070[56789]*.nc"
  files = systemfunc("ls " + histPath)

  if (any(ismissing(files))) then
    print ("    WARNING: No files found for ... " + histPath)
    delete(files)
    continue
  end if 
  
  nFiles = dimsizes(files)
  
  if (nFiles .le. 0) then
    print("    WARNING: No files found for ... " + histPath)
    delete(files)
    continue
  end if
  
  
  do f = 0, nFiles-1
    print("    " + files(f))
    history = addfile(files(f), "r")
    
    date = history->date
    nDates = dimsizes(date)
    
    slev = history->tpgpaltitude
    slev@_FillValue = -999
    slev = where(slev .eq. max(slev), slev@_FillValue, slev)
    tmp = dimsizes(slev)
    nsLevs = tmp(1)
    delete(tmp)
    
    slat = history->tplatitude
    slat@_FillValue = -999
    slat = where(slat .eq. max(slat), slat@_FillValue, slat)

    slon = history->tplongitude
    slon@_FillValue = -999
    slon = where(slon .eq. max(slon), slon@_FillValue, slon)
    
    print("      nDates :  " + nDates)
    print("      lev :  " + min(slev) + "  " + max(slev))
    print("      lat :  " + min(slat) + "  " + max(slat))
    print("      lon :  " + min(slon) + "  " + max(slon))
    
    ; Iterpolate to the verical grid
    sT = history->ktemp
    sT@_FillValue = -999
    sT = where(sT .eq. max(sT), sT@_FillValue, sT)
    
    do i = 0, nDates-1
      doy = floattoint(date(i) - (floor(date(i) / 1000) * 1000))
      
      if ((i % 1000) .eq. 0) then
        print("        " + i + "  " + doy)
      end if
      
      do k = 0, nsLevs - 1
        if (.not. ismissing(sT(i, k))) then
          plev = slev(i, k)
          plat = slat(i, k)
          plon = slon(i, k)
        
          ; The longitude lookup can't handle values between the largest lon and 360,
          ; so move them around.
          if (plon .gt. max(lon) + dlon/2) then
            plon = 0.
          end if

          if (plon .gt. max(lon)) then
            plon = doubletofloat(max(lon))
          end if
     
          nPoints({doy}, {plev}, {plat}, {plon}) = nPoints({doy}, {plev}, {plat}, {plon}) + 1
          T({doy}, {plev}, {plat}, {plon}) = T({doy}, {plev}, {plat}, {plon}) + sT(i, k)
        end if
      end do
    end do
    
    delete(slev)
    delete(slat)
    delete(slon)
    delete(sT)
    delete(date)
    delete(history)
  end do
    
  delete(files)
  
  ; Calculate the average temperature
  print("")
  print("  Calculating averages ...") 
  nPoints = where(nPoints .gt. 0, nPoints, nPoints@_FillValue)
  T = where(nPoints .gt. 0, T / nPoints, T@_FillValue)
  
  print("")
  print("  Writing Results ...")
  system("rm analysis/saber_T_wg.nc") 
  outfile               = addfile("analysis/saber_T_wg.nc", "c")
  outfile->gw           = gw
  outfile->nPoints      = nPoints
  outfile->T            = T

  print("")
  print(T&time + "  " + dim_avg(T(:, {82.}, {70.}, :)))
  print("")
  
  tmp = dim_avg_Wrap(dim_avg_Wrap(T(:, :, {60.:80.}, :)))
  print(T&lev + "  " + dim_avg(tmp(lev|:, time|:)))
  delete(tmp)

end


; Main
WACCM_RUN = getenv("WACCM_RUN")

if (ismissing(WACCM_RUN)) then
  WACCM_RUN = "results"
end if
print("Using WACCM_RUN: " + WACCM_RUN)
print("")


;gridSaber("analysis/"+WACCM_RUN+"/crice_alb_265_h3_06_cips.nc")
;gridSaber("analysis/old/pmc_rt_smin/crice_alb_265_h3_06_cips.nc")

;Note ptmp Samba share on Cadfael must be mounted on machine running this script.
;Input file just gives lat, lon, lev, and gw (Gaussian weights for latitude bands)
gridSaber("/Volumes/ptmp/mills/waccm/dust-sulfur/run/expuT125E01/expuT125E01.cam2.h3.1997-03-01-00000.nc")
