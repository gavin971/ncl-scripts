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
RHO_H2SO4 = 1.923                 ; density of pure H2SO4 (l) at 200K [ g / cm3 ]
RHO_DUST = 2.0                    ; density of meteoritic dust [ g / cm3 ]

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "ncl-scripts/getHistoryPath.ncl"

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

;=========================================

undef("setupbins")
procedure setupbins(nbin:integer, rmin:float, rmrat:float, rho:float, \
     rmass[*]:float, r[*]:float, dr[*]:float, rmassup[*]:float, dm[*]:float, \
     vol[*]:float, rup[*]:float, rlow[*]:float)

local cpi, i, vrfact

begin

  ; Allocate arrays for the return variables.
  r = new(nbin, float)
  r!0="bins"
  rmass = r
  dr = r
  
  ; NOTE: These aren't currently being returned, but are calculated
  ; in setupbins.
  rmassup = r
  dm = r
  vol = r
  rup = r
  rlow = r
  
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

;=========================================

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

;=========================================

; Get the number density for each bin.
undef("getN")
procedure getN(history:file, shortName:string, rmass:float, arho:float, N:float, skipbin:integer)
  
;  rmass(nBins) = mass of 1 particle in each bin (g)
;  arho(:,:,:,:) = atmospheric density (g/cm3)
;  N( ,b) = number density of particles in each bin (#/cm3) 
;  skipbin = bin to skip and interpolate (-1 to not skip any)
local b, binName, nBins, mmr

begin
  
  ; Iterate over the bins
  nBins = dimsizes(rmass)
    
;  ndname=shortName+"ND"
;  nd=history->$ndname$
  
;  ndsum=nd
;  ndsum(:,:,:,:)=0.0
  
  do b = 0, nBins-1
    if (b.ne.skipbin) then
      binName = shortName + sprinti("%2.2i", b+1)

      ; Get the mass mixing ratio.
      mmr = history->$binName$

      ; Compute the particle concentration.
      N(:, :, :, :, b) = (mmr * arho) / rmass(b)
;      ndsum(:,:,:,:)=ndsum(:,:,:,:)+N(:, :, :, :, b)
      delete(mmr)
      print("bin "+b+" max="+max(N(:, :, :, :, b))+" min="+min(N(:, :, :, :, b)))
    end if
  end do

; HACK: interpolate missing bin 32: 
  if (skipbin.gt.0) then 
    N(:, :, :, :, skipbin) = sqrt(N(:, :, :, :, skipbin-1)*N(:, :, :, :, skipbin+1))
  end if

;  N(:, :, :, :, 31) = (nd-ndsum>0)
;  mmr = N(:, :, :, :, 31)*rmass(31)/arho
;  mmr@units="kg/kg"
;  mmr@long_name="mass in bin 32"
;  binName = shortName + "32"
;  history->$binName$=mmr
 
end 

;=========================================

undef("saveSize")
procedure saveSize(historyRun:string,prefix:string,skipbin:integer)

;  inputs: historyRun, i.e. "coag04"
;          prefix, i.e. "nopur2.cam2.h0.1995"
;          skipbin, missing bin to skip & interpolate (-1 to not skip)

local historyNames,rmin,cpi,rmassmin,writefile,hrc,nhrc,fc

begin

  historyPath = "/ptmp/mills/waccm/carma/run/" + historyRun +"/"
    
; Create array of files matching pattern
  historyNames = systemfunc("ls " + historyPath + " | grep " + prefix)
  print(historyNames)

  ; Setup the same bin structure used during the model run.
  nBins = 28
  rmin  = 2e-8 ; cm
  rmrat = 2.0
  rho   = RHO_DUST  ; g/cm-3
  
  rmass = new(nBins, float)
  rmass!0 = "bins"
  dustrad = rmass
  dustdr = rmass
  rmassup = rmass
  dm = rmass
  dustvol = rmass
  dustradup = rmass
  dustradlow = rmass
  
  print("")
  print("    Assuming ...")
  print("      nBins    =  " + nBins)
  print("      rmin     =  " + rmin)
  print("      rmrat    =  " + rmrat)
  print("      rho      =  " + rho)
  
  setupbins(nBins, rmin, rmrat, rho, rmass, dustrad, dustdr, rmassup, dm, dustvol, dustradup, dustradlow)
  
  rmass@long_name  = "bin center mass"
  rmassup@long_name = "upper bin boundary mass"
  dustrad@long_name = "dust bin mean (volume-weighted) radius"
  dustradup@long_name = "dust upper bin boundary radius"
  dustdr@long_name = "dust bin radius spacing"
  dustradlow@long_name = "dust lower bin boundary radius"
  dustvol@long_name = "dust bin mean volume"
  dm@long_name = "bin mass spacing"
 
  rmass@units = "g"
  rmassup@units = "g"
  dustrad@units = "cm"
  dustradup@units = "cm"
  dustdr@units = "cm"
  dustradlow@units = "cm"
  dustvol@units = "cm3"
  dm@long_name = "g"
    
  
  nFiles=dimsizes(historyNames)
  do ifile = 0, nFiles-1
    print("")
    print("file " + (ifile+1) + " of " + nFiles + " :")
    print("   " + (/ historyNames(ifile) /))
    
    
    ; copy history file:
;    allfilename=historyPath+historyNames(ifile)
    allfilename=historyPath+historyNames(ifile)
    print(allfilename)
    writefile=addfile(allfilename,"w")

    ; Find the last time.
    time = writefile->time
    nFTimes = dimsizes(time)
    lev = writefile->lev
    ilev = writefile->ilev
    nLevs = dimsizes(lev)
    lat = writefile->lat
    nLats = dimsizes(lat)
    lon = writefile->lon
    nLons = dimsizes(lon)

    q     = writefile->Q
    temp  = writefile->T
            
    if (.not.(isfilevar(writefile,"RMASS"))) then
      print("writing aerosol properties")
      writefile->RMASS=rmass
      writefile->RMASSUP=rmassup
      writefile->DUSTRAD=dustrad
      writefile->DUSTRADUP=dustradup
      writefile->DUSTDR=dustdr
      writefile->DUSTRADLOW=dustradlow
      writefile->DUSTVOL=dustvol
      writefile->DM=dm
    end if
    
    if (isfilevar(writefile,"PRESSURE")) then
      P=writefile->PRESSURE
      arho=writefile->RHOAIR
    else
    
    ; Compute the air mass and density.
    ;
    ; NOTE: Convert arho and amass are in cgs units.
      P               = q ; get metadata and array structure from humidity
      P@units         = "Pa"
      P@long_name     = "Pressure"

      arho            = q ; get metadata and array structure from humidity
      arho@units      = "g/cm3"
      arho@long_name  = "air density"

      print("calling setupair")
      setupair(writefile, P, arho) ; calculate pressure and density
      arho = arho(:,:,:,:) * 1e3 / 1e6
      
      writefile->RHOAIR=arho
      writefile->PRESSURE=P
      
    end if
    

    ; Get the Number density (#/cm-3) for each bin at the 
    ; specified time and location.
 
    if (isfilevar(writefile,"DUST01") .and. (.not.(isfilevar(writefile,"DUSTNDBIN")))) then
      print("Calculating the dust number density ...")
      dustndbin = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
      dustndbin!0="time"
      dustndbin&time=time
      dustndbin!1="lev"
      dustndbin&lev=lev
      dustndbin!2="lat"
      dustndbin&lat=lat
      dustndbin!3="lon"
      dustndbin&lon=lon
      dustndbin!4="bins"   
      dustndbin@long_name="dust number density per bin"
      dustndbin@units="/cm3"
      getN(writefile, "DUST", rmass, arho, dustndbin, skipbin)
      writefile->DUSTNDBIN=dustndbin
    end if
  
            
    delete(time)
  end do
end

; Main

saveSize("125_lt_gw2_28","125_lt_gw2_28.cam2",-1) 

