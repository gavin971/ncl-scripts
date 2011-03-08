; Nuclear Soot

load "../core.ncl"
load "soot-ic-mass.ncl"

; Generate the initial tracer distribution.
;
; NOTE: The values for nbin, rmin, rmrat. and rho should match those
; in CARMA (aerrad.h and defineaer.f).
undef("newicBCIndPak")
procedure newicBCIndPak(srcIcPath:string, srcIcFile:string, icFileName:string, mass_tg:double)

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
;  mmr = 1.744568e-05 * (mass_tg / 5.0) * (3.8*5.0) / (50 * 1.9 * 2.5)
  mmr = 1.748958127e-05 * (mass_tg / 5.0) * (3.8*5.0) / (50 * 1.9 * 2.5)
  
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
  ;   input in two triangular regions over Pakistan and India
  ;   Apex of Pakistan triangle at 72.5E, 33.16N
  ;   300-150 mbar
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
  
  totalMass = 0.0d0 
  calculateMass(icPath, name, totalMass)
  
  mmr = mmr * (mass_tg / totalMass)
  print("new mmr="+mmr)
  
  do i =0,4
    tmp(0, {300:150}, j1-i, l1-i:l1+i) = mmr
    tmp(0, {300:150}, j2+i, l2-i:l2+i) = mmr
  end do
  
  ; Write the field.
  icFile->$name$ = (/ tmp /)
  
;
; For a nicer looking netCDF, create a "new line" character.
; This is not necessary.
;
  nl = integertochar(10)  ; newline character
;
; Define global attributes.
;
; globalAtt can be of any type. Here logical is used by convention.
;
  globalAtt             = True

  globalAtt@creation_date     = nl+\
      systemfunc("date")
      
  globalAtt@creator = nl+\
      "Mike Mills, mmills@ucar.edu"
  globalAtt@title = nl+\
      "IC file for India-Pakistan nuclear war study with 5Tg black carbon on May 15, 2012" ;
  globalAtt@history = nl+\
       "Adapted from "+icFileName+nl+\
       "SOOT01 field added by programs soot-ic.ncl and soot.ncl"
  globalAtt@notes = nl+\
       "5 Tg of black carbon input between 300 and 150 mbar in two triangular regions"+nl+\ 
       "over Pakistan and India. Apex of Pakistan triangle at 72.5E, 33.16N."+nl+\
       "SOOT01 supplied with a constant mass mixing ratio of "+mmr ;

  fileattdef( icFile, globalAtt )

end

