; Advection Test Case #1
; CAM-offline 3.1.3
;

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


undef("closest_lon")
function closest_lon(lon:numeric, lons[*])

; Longitudes wrap around, so a little care must be taken in the 360 -> 0 area.
;
; Parameters
;   (IN) lons
;   A vector of longitudes.
;
;   (IN) lon
;    The longitude to be searched for in the lons.

local iLon, nLons

begin

  iLon = closest_val(lon, lons)
  nLons = dimsizes(lons)

  ; If the value is great than the last entry, closest_val returns the
  ; number of dimensions.
  if (iLon .eq. nLons) then
    if ((min(lons) + 360. - lon) .gt. (lon - max(lons))) then
      iLon = 0
    else
      iLon = nLons - 1
    end if
  end if
  
  return(iLon)
end


undef("column_mmr")
function column_mmr(k0:integer, histTT[*][*]:numeric, histdp[*][*]:numeric)
  
local column, dims, k, nLevs, nPoints

begin
  ; Determine the dimensions.
  dims = dimsizes(histTT)
  nLevs = dims(0)
  nPoints = dims(1)
  
  ; Create an array to return.
  column = new(nPoints, double)
  column = 0.
  
  ; We get vertical motion that isn't wanted, so sum over the whole column
  ; and then assume all the mass is at ~k0.
  do k = 0, nLevs - 1
    column(:) = column(:) + histTT(k, :) * histdp(k, :)
  end do
    
  ; Express as a mass mixing ratio at k0.
  column = column / histdp(k0, :)
  
  return(column)
end


undef("lat_slice")
function lat_slice(down[*]:numeric, up[*]:numeric)
  
local slice, nLats, nPoints

begin

  ; Determine the dimensions.
  nLats = dimsizes(down)
  nPoints = 2 * nLats - 2
  
  ; Create an array to return.
  slice = new(nPoints, double)
  
  ; Go from north pole to south pole and then back to the north pole.
  ; This requires flipping the down array.
  slice(0: nLats-1) = down(::-1)
  slice(nLats:nPoints-1) = up(1:nLats-2)

  ; We need to create single values for the pole points.
  slice(0) = (up(nLats-1) + down(nLats-1)) / 2.
  slice(nLats-1) = (up(0) + down(0)) / 2.
  
  return(slice)
end



undef("createmetfile")
function createmetfile(srcIc:file, metPath:string, nDays:integer, dTime:integer)

; Creates a meteorological (met) data file based upon the data that is in the
; initial conditions (ic) file. The grid spacing (lev, lat and lon) will be
; taken from the ic file. Additional records will be added based upon ndays
; and dTime, but the fields exceot date and datesec will all have the same
; value as the ic file.
;
; The following dimensions will be transfered from the ic file to the met file:
;   time
;   lev
;   lat
;   lon
;
; The following fields will be transfered from the ic file to the met file:
;   lev(lev)
;   lat(lat)
;   lon(lon)
;   date(time)
;   datesec(time)
;   PS(time, lat, lon)
;   TAUX(time, lat, lon)
;   TAUY(time, lat, lon)
;   QFLX(time, lat, lon)
;   SHFLX(time, lat, lon)
;   T(time, lev, lat, lon)
;   U(time, lev, lat, lon)
;   V(time, lev, lat, lon)
;
; Parameters
;   (IN) srcIc
;   The NETCDF file containing the initial conditions from with this file will
;   be based.
;
;   (IN) metPath
;   The full path to the new file to be created.
;
;   (IN) nDays
;   The number of days for which to generate the met data.
;
;   (IN) dTime
;   The periodicity of the met data (seconds).
;
; Returns
;   A pointer to the newly created file.

local secsPerDay, metPath, nPerDay, nTime, dstMet, year, month, day, \
  startJDay, date, datesec, gday, iday, jday, i, j, vars, QFLX, SHFLX, TAUX, \
  TAUY, tmp, tmp2, startDate, startTime, U, V, US, VS

begin

  ; Constant(s)
  secsPerDay    = 24*60*60

  ; Remove any preexisting file.
  print("  Deleting " + metPath)
  system("/bin/rm " + metPath)

  ; Create the file.
  print("  Creating " + metPath)
  dstMet = addfile(metPath, "c")
  
  ; Create an unlimited time dimension
  filedimdef(dstMet, "time", -1, True)
  
  ; Determine the the data and time fields ...
  startDate = (/ srcIc->date(0) /)
  startTime = (/ srcIc->time(0) /)
  nPerDay = secsPerDay / dTime
  nTime = nDays * nPerDay
  
  ; We need to work with julian dates.
  print("    Creating " + nTime + " records starting on " + startDate)
  year = startDate / 10000
  month = (startDate % 10000) / 100
  day = startDate % 100
  jday = greg2jul(year, month, day, 0)

  ; Create the entries for the date and datesec variables.
  print("    Adding time, date, datesec ...")
  time = new(nTime, typeof(srcIc->time))
  date = new(nTime, "integer")
  datesec = new(nTime, "integer")

  do i = 0, nDays - 1
    gday = jul2greg(jday + i)
    iday = gday(0) * 10000 + gday(1) * 100 + gday(2)
    
    do j = 0, nPerDay - 1
      time(i * nPerDay + j) = startTime + i + (j * dTime) / (int2flt(secsPerDay))
      date(i * nPerDay + j) = iday
      datesec(i * nPerDay + j) = j * dTime
    end do
  end do

  ; Write the new time, date and datesec fields.
  filevardef(dstMet, "time", typeof(srcIc->time), (/ "time" /))
  filevarattdef(dstMet, "time", srcIc->time)
  dstMet->time = (/ time /)

  filevardef(dstMet, "date", typeof(srcIc->date), (/ "time" /))
  filevarattdef(dstMet, "date", srcIc->date)
  dstMet->date = (/ date /)
  
  filevardef(dstMet, "datesec", typeof(srcIc->datesec), (/ "time" /))
  filevarattdef(dstMet, "datesec", srcIc->datesec)
  dstMet->datesec = (/ datesec /)

  ; Variables that are independent of time can be copied over directly.
  print("    Adding lev, lat, lon, hyam, hybm, p0, slat, slon ...")
  dstMet->lev 	= srcIc->lev
  dstMet->lat 	= srcIc->lat
  dstMet->lon 	= srcIc->lon
  dstMet->hyam 	= srcIc->hyam
  dstMet->hybm 	= srcIc->hybm
  dstMet->P0 		= srcIc->P0
  dstMet->slat 	= srcIc->slat
  dstMet->slon 	= srcIc->slon

  ; We want to copy the following data, but increase the number of dimensions
  ; for time to the number of times needed for nDays and dTime.
  
  ; These are the surface variable that depend on time, lat and lon.
  print("    Adding PS, PHIS ...")
  vars = (/ "PS", "PHIS" /)
  do i = 0, dimsizes(vars) - 1
    filevardef(dstMet, vars(i), typeof(srcIc->$vars(i)$), (/"time","lat","lon"/))
    filevarattdef(dstMet, vars(i), srcIc->$vars(i)$)

    tmp = new(dimsizes(dstMet->$vars(i)$), typeof(dstMet->$vars(i)$))
    tmp2 = (/ srcIc->$vars(i)$(0, :, :) /)
    do j = 0, nTime - 1
      tmp(j, :, :) = (/ tmp2 /)
    end do
    dstMet->$vars(i)$ = (/ tmp /)
    delete(tmp)
    delete(tmp2)
  end do
  
  ; For some reason, QFLX, SHFLX, TAUX, and TAUY are not in the initial
  ; conditions file, so create them here.
  ;
  ; NOTE: It is really slow to write directly to the NETCDF file, so initialize
  ; the variables first and then write them out.
  print("    Adding QFLX, SHFLX, TAUX, TAUY ...")
  filevardef(dstMet, "QFLX", "double", (/"time","lat","lon"/))
  dstMet->QFLX@units        = "kg/kg"
  dstMet->QFLX@long_name    = "Surface water flux"
  QFLX = new(dimsizes(dstMet->QFLX), typeof(dstMet->QFLX))
  QFLX = 0.
  dstMet->QFLX					= (/ QFLX /)
  delete(QFLX)

  filevardef(dstMet, "SHFLX", "double", (/"time","lat","lon"/))
  dstMet->SHFLX@units        = "W/m2"
  dstMet->SHFLX@long_name    = "Surface sensible heat flux"
  SHFLX = new(dimsizes(dstMet->SHFLX), typeof(dstMet->SHFLX))
  SHFLX = 0.
  dstMet->SHFLX					= (/ SHFLX /)
  delete(SHFLX)
 
  filevardef(dstMet, "TAUX", "double", (/"time","lat","lon"/))
  dstMet->TAUX@units        = "N/m2"
  dstMet->TAUX@long_name    = "zonal surface stress"
  TAUX = new(dimsizes(dstMet->TAUX), typeof(dstMet->TAUX))
  TAUX = 0.
  dstMet->TAUX					= (/ TAUX /)
  delete(TAUX)
  
  filevardef(dstMet, "TAUY", "double", (/"time","lat","lon"/))
  dstMet->TAUY@units        = "N/m2"
  dstMet->TAUY@long_name    = "meridional surface stress"
  TAUY = new(dimsizes(dstMet->TAUY), typeof(dstMet->TAUY))
  TAUY = 0.
  dstMet->TAUY					= (/ TAUY /)
  delete(TAUY)
  
  ; These variables depend on time, lev, lat and lon.
  ;
  ; T is odd in that the fields are in a different order in the initial
  ; conditions file than they are in the meteorolgy file.
  print("    Adding Q, T ...")
  delete(vars)
  vars = (/ "Q", "T" /)
  do i = 0, dimsizes(vars) - 1
    filevardef(dstMet, vars(i), typeof(srcIc->$vars(i)$), (/"time","lev","lat","lon"/))
    filevarattdef(dstMet, vars(i), srcIc->$vars(i)$)

    tmp = new(dimsizes(dstMet->$vars(i)$), typeof(dstMet->$vars(i)$))
    tmp2 = (/ srcIc->$vars(i)$(time|0, lev|:, lat|:, lon|:) /)
    do j = 0, nTime - 1
      tmp(j, :, :, :) = (/ tmp2 /)
    end do
    dstMet->$vars(i)$ = (/ tmp /)
    delete(tmp)
    delete(tmp2)
  end do
  delete(vars)

  print("    Adding U, V, US, VS ...")
  filevardef(dstMet, "U", "double", (/"time","lev","lat","lon"/))
  dstMet->U@units        = "m/s"
  dstMet->U@long_name    = "Zonal wind"
  U = new(dimsizes(dstMet->U), typeof(dstMet->U))
  U = 0.
  dstMet->U					= (/ U /)
  delete(U)
  
  filevardef(dstMet, "V", "double", (/"time","lev","lat","lon"/))
  dstMet->V@units        = "m/s"
  dstMet->V@long_name    = "Meridional wind"
  V = new(dimsizes(dstMet->V), typeof(dstMet->V))
  V = 0.
  dstMet->V					= (/ V /)
  delete(V)
  
  filevardef(dstMet, "US", "double", (/"time","lev","slat","slon"/))
  dstMet->US@units        = "m/s"
  dstMet->US@long_name    = "Zonal wind, staggered"
  US = new(dimsizes(dstMet->US), typeof(dstMet->US))
  US = 0.
  dstMet->US				= (/ US /)
  delete(US)
  
  filevardef(dstMet, "VS", "double", (/"time","lev","lat","slon"/))
  dstMet->VS@units        = "m/s"
  dstMet->VS@long_name    = "Meridional wind, staggered"
  VS = new(dimsizes(dstMet->VS), typeof(dstMet->VS))
  VS = 0.
  dstMet->VS				= (/ VS /)
  delete(VS)
  
  return(dstMet)
end


undef("defaultatmosphere")
procedure defaultatmosphere(dstMet:file)

; Sets up a default atmosphere in the met file. The atmosphere will have the
; folowing characteristics:
;
;   - isothermal
;   - constant surface pressure
;   - constant geopotential height
;   - no surface stress
;   - no sensible heat flux
;   - exponetial water vapor concentration
;   - no water vapor flux
;
; Parameters
;   (IN) dstMet
;   The met file to modify.

local dQ, lev, log_minQ, minQ, maxQ, Q, QFLX, SHFLX, TAUX, TAUY, T

begin

  ; Read in some basic fields: lev, lat, lon, slat, slon
  lev = dstMet->lev
  
	print("  Setting Default Atmosphere ...")

	print("    Setting PS ...")
	PS = new(dimsizes(dstMet->PS), typeof(dstMet->PS))
	PS = 100000.
	dstMet->PS					= (/ PS /)
	delete(PS)

	print("    Setting PHIS ...")
	PHIS = new(dimsizes(dstMet->PHIS), typeof(dstMet->PHIS))
	PHIS = 0.
	dstMet->PHIS					= (/ PHIS /)
	delete(PHIS)

	; NOTE: Setting Q to a constant 1e-6 caused the run to fail during the
	; first step with NAN on some energy calculations. Try a gradient from
	; the top to the ground.
	print("    Setting Q ...")
	minQ = 1e-6
	maxQ = 1e-3
	log_minQ = log(minQ)
	nLevs = dimsizes(lev)
	dQ = (log(maxQ) - log_minQ) / (nLevs - 1)
	
	Q = new(dimsizes(dstMet->Q), typeof(dstMet->Q))

	do i = 0, nLevs - 1
		Q(:, i, :, :)       = (/ exp(log_minQ + i * dQ) /)
	end do
	
	dstMet->Q					= (/ Q /)
	delete(Q)

	print("    Setting QFLX ...")
	QFLX = new(dimsizes(dstMet->QFLX), typeof(dstMet->QFLX))
	QFLX = 0.
	dstMet->QFLX					= (/ QFLX /)
	delete(QFLX)
	
	print("    Setting SHFLX ...")
	SHFLX = new(dimsizes(dstMet->SHFLX), typeof(dstMet->SHFLX))
	SHFLX = 0.
	dstMet->SHFLX					= (/ SHFLX /)
	delete(SHFLX)
	
	print("    Setting TAUX ...")
	TAUX = new(dimsizes(dstMet->TAUX), typeof(dstMet->TAUX))
	TAUX = 0.
	dstMet->TAUX					= (/ TAUX /)
	delete(TAUX)

	print("    Setting TAUX ...")
	TAUY = new(dimsizes(dstMet->TAUY), typeof(dstMet->TAUY))
	TAUY = 0.
	dstMet->TAUY					= (/ TAUY /)
	delete(TAUY)
    
  print("    Setting T ...")
  T = new(dimsizes(dstMet->T), typeof(dstMet->T))
  T = 250.
  dstMet->T							= (/ T /)
  delete(T)
end


undef("solidybody")
procedure solidbody(dstMet:file, period:float, beta:float)

; Takes the met file and sets it up for solid body rotation. 
;
; Parameters
;   (IN) dstMet
;   The met file to modify.
;
;   (IN) period
;   The number of days for one rotation.
;
;   (IN) beta
;   The angle between the axis of rotation and the axis between the poles.

local pi, re, secsPerDay, srcIc, startDate, dstMet, time, lev, lat, lon, \
  Umax, QFLX, SHFLX, TAUX, TAUY, T, U, V, slat, slon, US, VS

begin

  ; Constant(s)
  pi      = 4.0 * atan(1.0)
  dg2rad  = pi / 180.0          ; rad deg^-1
  re      = 6.37122e06          ; m
  secsPerDay    = 24*60*60

	print("  Setting Solid Body Rotation ...")
	
  ; Read in some basic fields: lev, lat, lon, slat, slon
  lev = dstMet->lev
  lat = dstMet->lat
  lon = dstMet->lon
  slat = dstMet->slat
  slon = dstMet->slon
      
  ; Pick a speed that will get around the prime meridian in 30 days.
  VELmax = (2. * pi * re) / (period * secsPerDay)
    
  ; The equations for solid body ortation at an arbitary angle (beta) are:
  ;
  ;   u = U[cos(beta)cos(lat) + sin(beta)sin(lat)cos(lon)
  ;   v = -Usin(beta)sin(lon)
  ;
  ; From Williamson&Rasch[1989], p. 118
    
  dstMet@case_name = "solid body rotation, T=" + period + ", V=" + beta + ", Vmax=" + VELmax 
    
  print("    Setting U, US, V, VS ...")
  U  = new(dimsizes(dstMet->U),  typeof(dstMet->U))
  US = new(dimsizes(dstMet->US), typeof(dstMet->US))
  V  = new(dimsizes(dstMet->V),  typeof(dstMet->V))
  VS = new(dimsizes(dstMet->VS), typeof(dstMet->VS))
    
  sinbeta = sin(beta * dg2rad)
  cosbeta = cos(beta * dg2rad)
    
  do i = 0, dimsizes(lat) - 1
    sinlat = sin(lat(i) * dg2rad)
    coslat = cos(lat(i) * dg2rad)
      
    do j = 0, dimsizes(lon) - 1
    U(:, :, i, j)       = (/   VELmax * (cosbeta * coslat + sinbeta * sinlat * cos(lon(j) * dg2rad)) /)
    V(:, :, i, j)       = (/ - VELmax * sinbeta * sin(lon(j) * dg2rad) /)
    end do
    
    do j = 0, dimsizes(slon) - 1
      VS(:, :, i, j)       = (/ - VELmax * sinbeta * sin(slon(j) * dg2rad) /)
    end do
  end do
  
  dstMet->U             = (/ U /)
  delete(U)
  
  dstMet->V							= (/ V /)
  delete(V)

  dstMet->VS						= (/ VS /)
  delete(VS)

  do i = 0, dimsizes(slat) - 1
    sinlat = sin(slat(i) * dg2rad)
    coslat = cos(slat(i) * dg2rad)
    
    do j = 0, dimsizes(lon) - 1
      US(:, :, i, j)       = (/   VELmax * (cosbeta * coslat + sinbeta * sinlat * cos(lon(j) * dg2rad)) /)
    end do
  end do
  
  dstMet->US            = (/ US /)
  delete(US)
end


undef("syncic")
procedure syncic(icPath:string, metPath:string)

; Synchronizes the initial condition file with the meteorology file. This is
; done by setting the fields in the initial conditions file with the first
; values for the corresponding fields in the meteorology value.
;
; The folowing fields in the initial conditions will be reset:
;   date(time)
;   datesec(time)
;   PS(time, lat, lon)
;   PHIS(time, lat, lon)
;   QFLX(time, lat, lon)
;   SHFLX(time, lat, lon)
;   T(time, lev, lat, lon)

begin

  ; Open both files.
  dstIc = addfile(icPath, "w")
  dstMet = addfile(metPath, "r")

  ; Set the fields in the dstIc file to be the same as the first element
  ; in the Met file
  ;
  ; NOTE: Do I also need to do time?
  print("  Synchronizing PS, PHIS, TS, TS1, TS2, TS3, TS4")
  dstIc->PS(0, :, :)                 = (/ dstMet->PS(0, :, :) /)
  dstIc->PHIS(0, :, :)               = (/ dstMet->PHIS(0, :, :) /)
  dstIc->TS(0, :, :)                 = (/ dstMet->T(0, 0, :, :) /)
  dstIc->TS1(0, :, :)                = (/ dstMet->T(0, 0, :, :) /)
  dstIc->TS2(0, :, :)                = (/ dstMet->T(0, 0, :, :) /)
  dstIc->TS3(0, :, :)                = (/ dstMet->T(0, 0, :, :) /)
  dstIc->TS4(0, :, :)                = (/ dstMet->T(0, 0, :, :) /)
  ;dstIc->QFLX(0, :, :)              = (/ dst Met->QFLX(0, :, :) /)
  ;dstIc->SHFLX(0, :, :)             = (/ dstMet->SHFLX(0, :, :) /)
  ;dstIc->TAUX(0, :, :)              = (/ dstMet->TAUX(0, :, :) /)
  ;dstIc->TAUY(0, :, :)              = (/ dstMet->TAUY(0, :, :) /)

  print("  Synchronizing Q, T")
  dstIc->Q(0, :, :, :)              = (/ dstMet->Q(time|0, lat|:, lev|:, lon|:) /)
  dstIc->T(0, :, :, :)              = (/ dstMet->T(time|0, lat|:, lev|:, lon|:) /)
end


undef("calcmass")
function calcmass(histPath:string, tracerName:string)

local pi, dg2grad, g, re, histNames, hists, lev, lat, lon, nLevs, nLats, nLons, \
  tracerMMR, cellMass, columnMass, dlat, dlon, dx, dy, dydx, wgt, totalMass

begin

  ; Some constants we will need:
  pi      = 4.0 * atan(1.0)
  dg2rad  = pi / 180.0   ; rad deg^-1
  g       = 9.806               ; m s^-2 
  re      = 6.37122e06          ; m
  
  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, so we can do
  ; a time evolution of the mass.
  histNames = systemfunc("ls " + histPath)
  print("Reading files ...")
  print((/ histNames /))
  hists = addfiles(histNames, "r")
  ListSetType(hists, "cat")
  
  ; Read in the fields we will need.
  pdeldry = addfiles_GetVar(hists, histNames, "PDELDRY")
  tracerMMR = addfiles_GetVar(hists, histNames, tracerName)
  lev = tracerMMR&lev
  lat = tracerMMR&lat
  lon = tracerMMR&lon
  
  ; Determine some of the dimensions.
  nLevs = dimsizes(lev) 
  nLats = dimsizes(lat) 
  nLons = dimsizes(lon)
  
  ; For each grid point calculate the mass of tracer in the column.
  ;
  ; The amount of TRACER is specifed as a mass mixing ratio (kg/kg) to dry air.
  ; Assuming hydrostatic balance, the column mass at a pressure level is:
  ;    M = P / g
  ; The dry (i.e. no water vapor) pressure drop across the grid point is given
  ; by the field PDELDRY, so the mass of the column is the sum over all the
  ; levels of:
  ;    MASScolumn = SUMlev((PDELDRY / g) * TRACER)
  ;
  ; NOTE: By default, PDELDRY is output to the history files as an average value,
  ; but to check mass conservation, you want to work with instantaneous values.
  print("Calculating cell mass ...")
  cellMass = (pdeldry / g) * tracerMMR
  cellMass!0 = "time"
  cellMass!1 = "lev"
  cellMass&lev = lev
  cellMass!2 = "lat"
  cellMass&lat = lat
  cellMass!3 = "lon"
  cellMass&lon = lon

  print("Calculating column mass ...")
  columnMass = dim_sum_Wrap(cellMass(time|:, lat|:, lon|:, lev|:))

  ; Calculate the total mass of the tracer.
  ;
  ; The tracer mass at a surface grid point is the column mass multipled by
  ; the surface area of the grid point:
  ;    MASSgrid = GRIDarea * MASScolumn
  ; The total mass can then be computed adding up the mass of trace from each
  ; grid point:
  ;    MASStotal = SUMlat,lon(MASSgrid)
  
  ; NOTE: Use the wgt_areasum2 function to calculate the sum. The example2 at 
  ;   http://www.ncl.ucar.edu/Document/Functions/Built-in/wgt_areasum2.shtml
  ; shows how to do this.
  
  ; Determine a weighting factor for the sum proportional to the area of the
  ; grid point.
  
  print("Calculating total mass ...")
  dlon = fabs(lon(1) - lon(0))
  dlat = fabs(lat(1) - lat(0))
  
  ; We need to do some special handling to calculate these values to
  ; compensate for the poles and the equator.
  dlat2 = dlat / 2.
  dh = new((/ nLats /), typeof(lat))

  do i = 0, nLats - 1
    
    ; The normal handling, +- dlat / 2 in width.
    sdlat = sin((lat(i) + dlat2) * dg2rad) - sin((lat(i) - dlat2) * dg2rad)
    
    ; At the poles, we only want to dlat lat in one direction.
    if (lat(i) .eq. -90.)  then
      sdlat = 1 + sin((lat(i) + dlat2) * dg2rad)
    end if

    if (lat(i) .eq. 90.) then
      sdlat =  1 - sin((lat(i) - dlat2) * dg2rad)
    end if
    
    ; At the equator, the peak value is at the equator, so 2 times one side.
    if (lat(i) .eq. 0.)
      sdlat = 2. * (sin((lat(i) + dlat2) * dg2rad))
    end if

    dh(i) = re * sdlat
  end do
  
  dx = ((2. * pi) / nLons) * re
  dhdx = dh * dx
  
  print(dhdx)

  wgt = new((/ nLats, nLons /), double)
  wgt = conform(wgt, dhdx, 0)
  
  totalMass = wgt_areasum2(columnMass, wgt, 0)
  totalMass!0 = "time"
  tt = ispan(1, dimsizes(totalMass), 1)
  tt@long_name = "Time"
  tt@units = "timestep"
  totalMass&time = tt
  totalMass@units = "kg"
  totalMass@long_name = "Mass"
  
  print("intial mass  =  " + totalMass(0) + " (kg)")
  print("final mass   =  " + totalMass(dimsizes(totalMass) - 1) + " (kg)")
  print("delta mass   =  " + (totalMass(0) - totalMass(dimsizes(totalMass) - 1)) + " (kg)")
  print("min mass     =  " + min(totalMass) + " (kg)")
  print("max mass     =  " + max(totalMass) + " (kg)")
  
  one = new((/ 1, nLats, nLons /), typeof(tracerMMR))
  one = 1.0
  test = wgt_areasum2(one, wgt, 0)
  
  print("calced surface area :  " + 4. * pi * (re ^ 2))
  print("summed surface area :  " + test(0))
  
  return(totalMass)
end

function calculateMass(historyPath:string, tracerNames[*]:string)

local pi, dg2grad, g, re, historyNames, hists, lev, lat, lon, nLevs, nLats, nLons, \
  tracerMMR, cellMass, columnMass, dlat, dlon, dx, dy, dydx, wgt, mass, i, j, t, mass

begin

  ; Some constants we will need:
  pi      = 4.0 * atan(1.0)
  dg2rad  = pi / 180.0   ; rad deg^-1
  g       = 9.806               ; m s^-2 
  re      = 6.37122e06          ; m

  nTracers    = dimsizes(tracerNames)
  lastTotal   = 0
  
  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, so we can do
  ; a time evolution of the mass.
  print("  Reading data ...")
  historyNames = systemfunc("ls " + historyPath)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  
  ; Read in the fields we will need.
  pdeldry = addfiles_GetVar(history, historyNames, "PDELDRY")

  ; Create the output variable.
  nTimes = dimsizes(history[:]->time)
  mass = new((/ nTracers, nTimes /), "float")
  mass!0 = "tracer"
  mass&tracer = tracerNames
  mass!1 = "time"
  mass&time = history[:]&time
  mass@units = "kg"

  do t = 0, nTracers - 1
  
    ; Check for a special tracer name that totals all the other mass.
    print ("")
    
    if (tracerNames(t) .eq. "TOTAL") then

      print("  Calculating subtotal ...")
      mass(t, :) = mass(lastTotal, :)
      
      do t2 = lastTotal + 1, t-1
        mass(t, :) = mass(t, :) + mass(t2, :) 
      end do

    else
    
      print("  Calculating mass for " + tracerNames(t) + " ...")
      tracerMMR = addfiles_GetVar(history, historyNames, tracerNames(t))
  
      lev = tracerMMR&lev
      lat = tracerMMR&lat
      lon = tracerMMR&lon
    
      ; Determine some of the dimensions.
      nLevs = dimsizes(lev) 
      nLats = dimsizes(lat) 
      nLons = dimsizes(lon)
    
      ; For each grid point calculate the mass of tracer in the column.
      ;
      ; The amount of TRACER is specifed as a mass mixing ratio (kg/kg) to dry air.
      ; Assuming hydrostatic balance, the column mass at a pressure level is:
      ;    M = P / g
      ; The dry (i.e. no water vapor) pressure drop across the grid point is given
      ; by the field PDELDRY, so the mass of the column is the sum over all the
      ; levels of:
      ;    MASScolumn = SUMlev((PDELDRY / g) * TRACER)
      ;
      ; NOTE: By default, PDELDRY is output to the history files as an average value,
      ; but to check mass conservation, you want to work with instantaneous values.
      print("    cell mass ...")
      cellMass = (pdeldry / g) * tracerMMR
      cellMass!0 = "time"
      cellMass!1 = "lev"
      cellMass&lev = lev
      cellMass!2 = "lat"
      cellMass&lat = lat
      cellMass!3 = "lon"
      cellMass&lon = lon
    
      print("    column mass ...")
      columnMass = dim_sum_Wrap(cellMass(time|:, lat|:, lon|:, lev|:))
    
      ; Calculate the total mass of the tracer.
      ;
      ; The tracer mass at a surface grid point is the column mass multipled by
      ; the surface area of the grid point:
      ;    MASSgrid = GRIDarea * MASScolumn
      ; The total mass can then be computed adding up the mass of trace from each
      ; grid point:
      ;    MASStotal = SUMlat,lon(MASSgrid)
      
      ; NOTE: Use the wgt_areasum2 function to calculate the sum. The example2 at 
      ;   http://www.ncl.ucar.edu/Document/Functions/Built-in/wgt_areasum2.shtml
      ; shows how to do this.
      
      ; Determine a weighting factor for the sum proportional to the area of the
      ; grid point.
      
      print("    total mass ...")
      dlon = fabs(lon(1) - lon(0))
      dlat = fabs(lat(1) - lat(0))
      
      ; We need to do some special handling to calculate these values to
      ; compensate for the poles and the equator.
      dlat2 = dlat / 2.
      dh = new((/ nLats /), typeof(lat))
    
      do i = 0, nLats - 1
        
        ; The normal handling, +- dlat / 2 in width.
        sdlat = sin((lat(i) + dlat2) * dg2rad) - sin((lat(i) - dlat2) * dg2rad)
        
        ; At the poles, we only want to dlat lat in one direction.
        if (lat(i) .eq. -90.)  then
          sdlat = 1 + sin((lat(i) + dlat2) * dg2rad)
        end if
    
        if (lat(i) .eq. 90.) then
          sdlat =  1 - sin((lat(i) - dlat2) * dg2rad)
        end if
        
        ; At the equator, the peak value is at the equator, so 2 times one side.
        if (lat(i) .eq. 0.)
          sdlat = 2. * (sin((lat(i) + dlat2) * dg2rad))
        end if
    
        dh(i) = re * sdlat
      end do
      
      dx = ((2. * pi) / nLons) * re
      dhdx = dh * dx
      
      wgt = new((/ nLats, nLons /), double)
      wgt = conform(wgt, dhdx, 0)
      
      mass(t, :) = wgt_areasum2(columnMass, wgt, 0)
    end if
    
  end do

;  one = new((/ 1, nLats, nLons /), typeof(tracerMMR))
;  one = 1.0
;  test = wgt_areasum2(one, wgt, 0)
;
;  print("")  
;  print("Summation check ...")  
;  print("  calced surface area :  " + 4. * pi * (re ^ 2))
;  print("  summed surface area :  " + test(0))
  
  return(mass)
end


procedure plotMass(basePath:string, historyPattern:string, runNames[*]:string, \
  tracerNames[*]:string)
  
local nRuns, nTracers, historyPath, history, nTimes, mass, wks, res, resP, i

begin

  nRuns = dimsizes(runNames)
  nTracers = dimsizes(tracerNames)
  
  ; To figure out the size of the array that we will need, open the first
  ; set f history files.
  historyPath = basePath + runNames(0) + "/" + historyPattern
  historyNames = systemfunc("ls " + historyPath)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  nTimes = dimsizes(history[:]&time) 

  ; Create arrays to store all the mass data.
  mass = new((/ nRuns, nTracers, nTimes /), "float")
  mass!0 = "run"
  mass&run = runNames
  mass!1 = "tracer"
  mass&tracer = tracerNames
  mass!2 = "time"
  mass&time = history[:]&time
  mass@units = "kg"
  
  dmass = new((/ nRuns, nTracers, nTimes /), "float")
  dmass!0 = "run"
  dmass&run = runNames
  dmass!1 = "tracer"
  dmass&tracer = tracerNames
  dmass!2 = "time"
  dmass&time = history[:]&time
  dmass@units = "kg"
  
  ; Calculate the tracer masses for all the runs.
  do i = 0, nRuns - 1
    print("")
    print("")
    print("Run :  " + runNames(i))
    historyPath = basePath + runNames(i) + "/" + historyPattern
    mass(i, :, :) = calculateMass(historyPath, tracerNames)

    print("")
    do j = 0, nTracers - 1
      dmass(i, j, :) = mass(i, j, :) - mass(i, j, 0)
      
      print("  " + tracerNames(j) + " :  " + mass(i, j, nTimes-1) + " - " + \
        mass(i, j, 0) + " = " + dmass(i, j, nTimes-1))
    end do
  end do
  
  ; Create a panel plot of the mass difference.
  print("")
  print(" Plotting ...")
  wks = gsn_open_wks("pdf", "dmass")
  plot = new(nRuns, graphic)
  
  res = True
;  res@gsnMaximize           = True
  res@vpHeightF             = 1.0 / nRuns              ; change aspect ratio
  res@vpWidthF              = 1.0
  res@gsnDraw = False
  res@gsnFrame = False
  
  res@tmYLFormat = "0*+^e"
  res@tmYLPrecision = 3
  
  do i = 0, nRuns - 1
    res@gsnLeftString = runNames(i)
    res@gsnRightString = "kg"
    plot(i) = gsn_csm_xy(wks, dmass&time, dmass(i, nTracers - 1, :) / mass(i, nTracers-1, 0), res)
  end do
  
  resP = True
  resP@txString = "Change in Total Mass"
  gsn_panel(wks, plot, (/ nRuns, 1 /), resP)
  delete(wks)
end
