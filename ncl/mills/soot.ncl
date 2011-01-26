; Nuclear Soot

load "core.ncl"

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
;    tmp(0, {300:150}, {30}, {70}) = 1.2e-5
;    tmp(0, {300:150}, {30}, {70}) = 1.8031e-5 ; 5 Tg soot
    tmp(0, {300:150}, {30}, {70}) = 3.6062e-6 ; 1 Tg soot
;    tmp(0, 52:56, 31, 15) = 1.2e-5
    
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
procedure newic(srcIcPath:string, srcIcFile:string, icFileName:string, mass_tg:float)

begin

  icPath    = srcIcPath+"/"+icFileName
  
  ; Copy the initial conditions file.
  print("  Creating " + icPath)
  print("    from " + srcIcPath + "/" + srcIcFile)
  system("cp " + srcIcPath + "/" + srcIcFile + " " + icPath)

  ; Add initial conditions for Sulfate CN from a
  ; lognormal distribution.
  ;
  ; NOTE: Only add these values in the tropopause region
  ; This will eliminate cloud formation in other areas.
  print("  Adding particles ... ")
  
  ; 12 ppmm x 4 deg lat x 5 deg long x 300 - 150 mbar = 5 Tg Soot
;  mmr = 1.8031e-5 ; 5 Tg soot @ 4 deg lat x 5 deg long x 300 - 150 mbar
;  mmr = 3.6062e-6 ; 1 Tg soot @ 4 deg lat x 5 deg long x 300 - 150 mbar
;  mmr = 1.8031e-5 * (mass_tg / 5.0) * (20.0 / (dlat * dlon))
;  mmr = 1.8031e-5 * (mass_tg / 5.0) * (3.8*5.0) / (50 * 1.9 * 2.5)
  mmr = 1.744568e-05 * (mass_tg / 5.0) * (3.8*5.0) / (50 * 1.9 * 2.5)
  
  print("mmr="+mmr)
      
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
  
  ; Name the bin.
  name = "SOOT01"

  ; Create the field if it does not exist.
  if (.not.(isfilevar(icFile,name))) then
    filevardef(icFile, name, "double", (/"time","lev","lat","lon"/))
    icFile->$name$@units        = "kg/kg"
    icFile->$name$@long_name    = "CARMA particle " + name
  end if
  
  ; Set the mmr.
  tmp = 0.0
  ; Add initial mass mixing ratio distribution for nuclear soot here
  ;   input at (70E, 300-150 mbar, 30N)
  j1=closest_val(33.16,lat)
  j2=j1-9
  l1=closest_val(72.5,lon)
  l2=l1+3
  
  do i =0,4
    tmp(0, {300:150}, j1-i, l1-i:l1+i) = mmr
    tmp(0, {300:150}, j2+i, l2-i:l2+i) = mmr
  end do
  
  ; Write the field.
  icFile->$name$ = (/ tmp /)

end

