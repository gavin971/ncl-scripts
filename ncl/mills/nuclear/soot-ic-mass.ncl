load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

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

undef("calculateMass")
procedure calculateMass(icPath:string, binName:string, totalMass:double)

begin

  history = addfile(icPath, "r")

  ; Find the last time.
  time = history->time
  nTimes = dimsizes(time)
  lev = history->lev
  nLevs = dimsizes(lev)
  lat = history->lat
  nLats = dimsizes(lat)
  lon = history->lon
  nLons = dimsizes(lon)

  ; Setup the same bin structure used during the model run.
  nBins=1
  rmin  = 1e-7 ; cm
  rmrat = 2.0
  rho   = 1.0  ; g/cm-3

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

  ; Compute the air mass and density.
  ;
  ; NOTE: Convert arho and amass are in cgs units.
  print("")
  print("  Calculating air pressure and density ...")
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
  print("") 
  print("  Calculating bin mass ...")
  q         = history->Q
  mass = new((/ nTimes, nLevs, nLats, nLons /), double)
  mass!0 = "time"
  mass&time = time
  mass!1 = "lev"
  mass&lev = lev
  mass!2 = "lat"
  mass&lat = lat
  mass!3 = "lon"
  mass&lon = lon

  ; NOTE: To keep the memory footprint smaller, calculate the total
  ; mass for each bin one at a time.
;  do b = 0, nBins-1
;    binName = "SOOT" + sprinti("%2.2i", bin+1)
    print("    " + binName)

    ; Get the mass mixing ratio.
    mmr = history->$binName$

    ; Compute the particle concentration.
    ; NOTE: Ignore the mass at the bottom level, since this is where the mass that
    ; sediments will collect. Perhaops we should change the model to get rid of this
    ; mass.
    mass             = (mmr(:, :, :, :) * (1 + q(:, :, :, :)) * (pdel(:, :, :, :) / G))        
    columnMass       = dim_sum_Wrap(mass(time|:, lat|:, lon|:, lev|0:nLevs-2))       
    totalMass  = wgt_areasum2(columnMass(:, :, :), wgt, 0)
    totalMass  = totalMass/1e9

    delete(mmr)
    delete(columnMass)
    print(totalMass)
;  end do

  delete(mass)
  delete(pdel)
  delete(q)           

;  print(totalMass)

end

;  MODEL_PATH = getenv("MODEL_PATH")

;  if (ismissing(MODEL_PATH)) then
;    MODEL_PATH  = "/waccm/mmills"
;  end if

;  icPath=MODEL_PATH+"/ccsm/inputdata/atm/waccm/ic/"+\
;         "b40.rcp4_5.2deg.wcm.carma.bc5tg.IndPak.cam2.i.2012-05-15-00000.nc"
icPath="/u/mjmills2/csm/b40.rcp4_5.2deg.wcm.carma.bc0tg.002/rest/2012-05-15-00000/"+\
       "b40.rcp4_5.2deg.wcm.carma.bc5tg.IndPak.002.cam2.i.2012-05-15-00000.nc"

totalMass=0.d0
bin=0
binName = "SOOT" + sprinti("%2.2i", bin+1)

calculateMass(icPath, binName, totalMass)
