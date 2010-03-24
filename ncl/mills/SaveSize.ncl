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
procedure setupbins(nbin:integer, rmin:float, rmrat:float, rho:float, rmass[*]:float, r[*]:float, dr[*]:float, rmassup[*]:float, dm[*]:float, vol[*]:float, rup[*]:float, rlow[*]:float)

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


; Compute sulfate mass density
undef("calcSulfDens")
procedure calcSulfDens(wtp:float,temp:float,sulfDens:float)
  
local DNWTP, DNC0, DNC1, ndim, iTime, iLev, iLat, iLon, den1, frac, i

begin

  DNWTP = (/0, 1, 5, 10, 20, 25, 30, 35, 40, 41, 45, 50, 53, 55, 56, 60, 65, 66, 70, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 100/)
  
  DNC0 = (/1, 1.13185, 1.17171, 1.22164, 1.3219, 1.37209, 1.42185, 1.4705, 1.51767, 1.52731, 1.56584, 1.61834, 1.65191, 1.6752, 1.68708, 1.7356, 1.7997, 1.81271, 1.86696, 1.89491, 1.9092, 1.92395, 1.93904, 1.95438, 1.98574, 2.00151, 2.01703, 2.03234, 2.04716, 2.06082, 2.07363, 2.08461, 2.09386, 2.10143, 2.10764, 2.11283, 2.11671, 2.11938, 2.12125, 2.1219, 2.12723, 2.12654, 2.12621, 2.12561, 2.12494, 2.12093/)
  
  DNC1= (/0, -0.000435022, -0.000479481, -0.000531558, -0.000622448, -0.000660866, -0.000693492, -0.000718251, -0.000732869, -0.000735755, -0.000744294, -0.000761493, -0.000774238, -0.00078392, -0.000788939, -0.00080946, -0.000839848, -0.000845825, -0.000874337, -0.000890074, -0.00089873, -0.000908778, -0.000920012, -0.000932184, -0.000959514, -0.000974043, -0.000988264, -0.00100258, -0.00101634, -0.00102762, -0.00103757, -0.00104337, -0.00104563, -0.00104458, -0.00104144, -0.00103719, -0.00103089, -0.00102262, -0.00101355, -0.00100249, -0.00100934, -0.000998299, -0.000990961, -0.000985845, -0.000984529, -0.000989315/)

  ndim=dimsizes(wtp)
  
  do iTime=0,ndim(0)-1
    print("iTime="+iTime)
    do iLev=0,ndim(1)-1
      do iLat=0,ndim(2)-1
        do iLon=0,ndim(3)-1
          i = 0
          do while (wtp(iTime,iLev,iLat,iLon) .gt. DNWTP(i))
            i=i+1
          end do

          sulfDens(iTime,iLev,iLat,iLon)=DNC0(i)+DNC1(i)*temp(iTime,iLev,iLat,iLon)

          if (i.ne.1 .and. wtp(iTime,iLev,iLat,iLon).ne.DNWTP(i)) then
            den1=DNC0(i-1)+DNC1(i-1)*temp(iTime,iLev,iLat,iLon)
            frac=(DNWTP(i)-wtp(iTime,iLev,iLat,iLon))/(DNWTP(i)-DNWTP(i-1))
            sulfDens(iTime,iLev,iLat,iLon) = den1 * frac + sulfDens(iTime,iLev,iLat,iLon) * (1.0-frac)
          end if
        end do
      end do
      print(" iLev="+iLev+": max="+max(sulfDens(iTime,iLev,:,:))+", min="+min(sulfDens(iTime,iLev,:,:)))
    end do
    print("iTime="+iTime+": max="+max(sulfDens(iTime,:,:,:))+", min="+min(sulfDens(iTime,:,:,:)))
  end do
  delete(DNWTP)
  delete(DNC0)
  delete(DNC1)
  delete(ndim)
  delete(iTime)
  delete(iLev)
  delete(iLat)
  delete(iLon)
  delete(den1)
  delete(frac)
  delete(i)
end


undef("calcSulfRad")
procedure calcSulfRad(history:file, mixsulnd:float, wtpct:float, sulfDens:float, coreName:string, rmass:float, arho:float, sulfRad:float, skipbin:integer)

local b, binName, nBins, mmr, coremass, cpi, voltot

begin
  cpi = 4./3. * PI

  dims=dimsizes(wtpct)
  ntime=dims(0)
  nlev=dims(1)
  nlat=dims(2)
  nlon=dims(3)
  
  ; Iterate over the bins
  nBins = dimsizes(rmass)
    
  do b = 0, nBins-1
;    print("bin "+b)
    if (b.ne.skipbin) then
      binName = coreName + sprinti("%2.2i", b+1)

      ; Get the mass mixing ratio.
      mmr = history->$binName$

      ; Compute core mass (g/particle)
      coremass = (mmr * arho) / (mixsulnd (:, :, :, :, b)>1e-30)

      voltot = coremass/RHO_DUST + (rmass(b)-coremass) * 100.0/wtpct / sulfDens

      sulfRad(:, :, :, :, b) = (voltot/cpi)^(1./3.)
      print("bin "+b+", sulfrad max="+max(sulfRad(:, :, :, :, b))+" min="+min(sulfRad(:, :, :, :, b)))
  
      delete(mmr)
    end if
    
    if (skipbin.gt.0) then 
      sulfRad(:, :, :, :, skipbin) = sqrt(sulfRad(:, :, :, :, skipbin-1)*sulfRad(:, :, :, :, skipbin+1))
    end if
    
  end do
  
end


undef("saveSulfDens")
procedure saveSulfDens(historyRun:string,year:string)

;  inputs: historyRun, i.e. "coag04"
;          year, i.e. "1995"

local historyNames,nFiles,ifile,history,writefile,wtpct,temp,sulfDens,hrc,nhrc,fc

begin
  historyPattern=historyRun+"/*.cam2.h0."+year+"*.nc"
; Create array of files matching pattern
  historyNames = systemfunc("ls " + historyPattern)    
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)
  
  nFiles=dimsizes(historyNames)
  do ifile = 0, nFiles-1
;    ifile = 17
    print("")
    print("    " + ifile + " :  " + (/ historyNames(ifile) /))

    history = addfile(historyNames(ifile), "r")

    wtpct = history->WTPCT
    temp  = history->T
    
    fc=stringtochar(historyNames(ifile))    
    
    writefile=addfile(historyRun+"/size."+chartostring(fc(nhrc:)),"c")   
;    writefile=addfile("sulfdens.nc","c")   
    
    sulfDens=wtpct
    sulfDens@long_name="sulfate aerosol mass density"
    sulfDens@units="g/cm3"
    print("      Calculating sulfate density ...")
    calcSulfDens(wtpct,temp,sulfDens)
    delete(wtpct)
    delete(temp)
    print("writing to sulfdens.nc ...")
;    writefile->SULFDENS=sulfDens   ; removed 2/21/08 (used before SULDEN saved from model)
    writefile->SULDEN=sulfDens
    print("deleting sulfDens")
    delete(sulfDens)
    delete(writefile)    
  end do
  delete(historyNames)
  delete(hrc)
  delete(nhrc)
end


undef("saveSize")
procedure saveSize(historyRun:string,year:string,skipbin:integer)

;  inputs: historyRun, i.e. "coag04"
;          year, i.e. "1995"
;          skipbin, missing bin to skip & interpolate (-1 to not skip)

local historyNames,rmin,cpi,rmassmin,readfile,writefile,hrc,nhrc,fc

begin
  historyPattern=historyRun+"/*.cam2.h0."+year+"*.nc"
; Create array of files matching pattern
  historyNames = systemfunc("ls " + historyPattern)
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)

  ; Setup the same bin structure used during the model run.
  nBins = 38
  rmin  = 2e-8 ; cm
  rmrat = 2.0
  rho   = RHO_DUST  ; g/cm-3
  
  rmass = new(nBins, float)
  rmass!0="bins"
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
  
  rmass@long_name  ="bin center mass"
  rmassup@long_name="upper bin boundary mass"
  dustrad@long_name="dust bin mean (volume-weighted) radius"
  dustradup@long_name="dust upper bin boundary radius"
  dustdr@long_name="dust bin radius spacing"
  dustradlow@long_name="dust lower bin boundary radius"
  dustvol@long_name="dust bin mean volume"
  dm@long_name="bin mass spacing"
 
  rmass@units  ="g"
  rmassup@units="g"
  dustrad@units="cm"
  dustradup@units="cm"
  dustdr@units="cm"
  dustradlow@units="cm"
  dustvol@units="cm3"
  dm@long_name="g"
    
  
  nFiles=dimsizes(historyNames)
  do ifile = 0, nFiles-1
    print("")
    print("    " + ifile + " :  " + (/ historyNames(ifile) /))

    history = addfile(historyNames(ifile), "r")

    ; Find the last time.
    time = history->time
    nFTimes = dimsizes(time)
    lev = history->lev
    ilev = history->ilev
    nLevs = dimsizes(lev)
    lat = history->lat
    nLats = dimsizes(lat)
    lon = history->lon
    nLons = dimsizes(lon)

    q     = history->Q
    wtpct = history->WTPCT
    temp  = history->T
    mixsulvmr 	= history->MIXSULMD
    pursulvmr 	= history->PURSULMD    
    soy      	= history->OCS    
    s 	       = history->S    
    so 	      = history->SO
    so2 	     = history->SO2
    so3 	     = history->SO3
    hso3      = history->HSO3   
    h2so4     = history->H2SO4    
    
    fc=stringtochar(historyNames(ifile))
    filen=historyRun+"/size."+chartostring(fc(nhrc:))
    print(filen)
    
    if (isfilepresent(filen)) then
      writefile=addfile(filen,"w")
    else
      writefile=addfile(filen,"c")
    end if
    
    if (isfilevar(history,"SULDEN")) then
      sulfDens=history->SULDEN    
    else
      if (isfilevar(writefile,"SULDEN")) then
        sulfDens=writefile->SULDEN
      else
        sulfDens=wtpct
        sulfDens@long_name="sulfate aerosol mass density"
        sulfDens@units="g/cm3"
        print("      Calculating sulfate density ...")
        calcSulfDens(wtpct,temp,sulfDens)
        writefile->SULDEN=sulfDens
      end if
    end if
    
    writefile->time=time
    writefile->lev=lev
    writefile->lat=lat
    writefile->lon=lon
    writefile->hyam 	= history->hyam
    writefile->hybm 	= history->hybm
    writefile->hyai 	= history->hyai
    writefile->hybi 	= history->hybi
    writefile->P0 	  = history->P0
    writefile->PS   	= history->PS

    writefile->RMASS=rmass
    writefile->RMASSUP=rmassup
    writefile->DUSTRAD=dustrad
    writefile->DUSTRADUP=dustradup
    writefile->DUSTDR=dustdr
    writefile->DUSTRADLOW=dustradlow
    writefile->DUSTVOL=dustvol
    writefile->DM=dm
    
;    sulfDens=wtpct
;    sulfDens@long_name="sulfate aerosol mass density"
;    sulfDens@units="g/cm3"
;    print("      Calculating sulfate density ...")
;    calcSulfDens(wtpct,temp,sulfDens)    
;    writefile->SULFDENS=sulfDens

    ; Compute the air mass and density.
    ;
    ; NOTE: Convert arho and amass are in cgs units.
    P               = q ; get metadata and array structure from humidity
    P@units         = "Pa"
    P@long_name     = "Pressure"

    arho            = q ; get metadata and array structure from humidity
    arho@units      = "g/cm3"
    arho@long_name  = "air density"

    setupair(history, P, arho) ; calculate pressure and density
    arho = arho(:,:,:,:) * 1e3 / 1e6

; (mol/mol)=(g/cm3)/(g/cm3)*(g/mol Air)/(g/mol H2SO4)
;   (28.964 g/mol Air)/(98.0734 g/mol H2SO4) =0.2953
    mixsulvmr=mixsulvmr(:,:,:,:)/arho(:,:,:,:)*0.2953
    mixsulvmr@units      = "mol/mol"
    mixsulvmr@long_name  = "pure sulfates equivalent vmr"
    
    pursulvmr=pursulvmr(:,:,:,:)/arho(:,:,:,:)*0.2953
    pursulvmr@units      = "mol/mol"
    pursulvmr@long_name  = "pure sulfates equivalent vmr"

    soy=soy+s+so+so2+so3+hso3+h2so4+mixsulvmr+pursulvmr
    soy@units      = "mol/mol"
    soy@long_name  = "total sulfur"

    ; Get the Number density (#/cm-3) for each bin at the 
    ; specified time and location.
    print("      Calculating the mixed sulfate number density ...")
    mixsulndbin = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
    mixsulndbin!0="time"
    mixsulndbin&time=time
    mixsulndbin!1="lev"
    mixsulndbin&lev=lev
    mixsulndbin!2="lat"
    mixsulndbin&lat=lat
    mixsulndbin!3="lon"
    mixsulndbin&lon=lon
    mixsulndbin!4="bins"   
    mixsulndbin@long_name="mixed sulfate number density per bin"
    mixsulndbin@units="/cm3"
    getN(history, "MIXSUL", rmass, arho, mixsulndbin, skipbin)
    writefile->MIXSULNDBIN=mixsulndbin
 
    print("      Calculating the dust number density ...")
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
    getN(history, "DUST", rmass, arho, dustndbin, skipbin)
    writefile->DUSTNDBIN=dustndbin
 
    sulfRad = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
    sulfRad!0="time"
    sulfRad&time=time
    sulfRad!1="lev"
    sulfRad&lev=lev
    sulfRad!2="lat"
    sulfRad&lat=lat
    sulfRad!3="lon"
    sulfRad&lon=lon
    sulfRad!4="bins"   
    sulfRad@long_name="mixed sulfate radius"
    sulfRad@units="cm"
    
    calcSulfRad(history, mixsulndbin, wtpct, sulfDens, "DUSSUL", rmass, arho, sulfRad, skipbin)
    writefile->SULFRAD=sulfRad
    writefile->RHOAIR=arho
    writefile->PRESSURE=P
    writefile->SOY=soy
    writefile->MIXSULVMR=mixsulvmr
    writefile->PURSULVMR=pursulvmr
    delete(time)
  end do

end

undef("addHybrid")
procedure addHybrid(historyRun:string,year:string)

;  inputs: historyRun, i.e. "coag04"
;          year, i.e. "1995"

local historyNames,rmin,cpi,rmassmin,readfile,writefile,hrc,nhrc,fc

begin
  historyPattern=historyRun+"/*.cam2.h0."+year+"*.nc"
; Create array of files matching pattern
  historyNames = systemfunc("ls " + historyPattern)
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)
  
  nFiles=dimsizes(historyNames)
  do ifile = 0, nFiles-1
    print("")
    print("    " + ifile + " :  " + (/ historyNames(ifile) /))

    history = addfile(historyNames(ifile), "r")
    fc=stringtochar(historyNames(ifile))
    filen=historyRun+"/size."+chartostring(fc(nhrc:))
    print(filen)
    writefile=addfile(filen,"w")   

    writefile->hyai  = history->hyai
    writefile->hybi 	= history->hybi   
  end do
end


; Main

;saveSulfDens("coag04","1995-01-01")
;addHybrid("coag04","1995-0[1-7]")
;addHybrid("coag04","1995-08-[01]")
;addHybrid("coag04","1995-08-24")

;saveSize("coag04","1995-08-29",31) ; skip bin 31 in coag04 directory
;saveSize("coag04","1995-09-23",31) ; skip bin 31 in coag04 directory

;saveSize("coag04-2","199",-1) ; include all bins in coag04-2 directory
;saveSize("coag04-2","1996-12-31",-1) ; include all bins in coag04-2 directory
;saveSize("coag04-2","1997",-1) ; include all bins in coag04-2 directory
;addRhoPressSOy("coag04-2","199")
saveSize("nopur2_ncar","1995",-1) ; include all bins in coag04-2 directory
