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
;  r = new(nbin, float)
;  r!0="bins"
;  rmass = r
;  dr = r
  
  ; NOTE: These aren't currently being returned, but are calculated
  ; in setupbins.
;  rmassup = r
;  dm = r
;  vol = r
;  rup = r
;  rlow = r
  
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
procedure getN(history:file, shortName:string, rmass:float, arho:float, N:float )
  
;  rmass(nBins) = mass of 1 particle in each bin (g)
;  arho(:,:,:,:) = atmospheric density (g/cm3)
;  N( ,b) = number density of particles in each bin (#/cm3) 
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
;      ndsum(:,:,:,:)=ndsum(:,:,:,:)+N(:, :, :, :, b)
      delete(mmr)
      print(binName+" max N="+max(N(:, :, :, :, b))+" min N="+min(N(:, :, :, :, b)))
  end do
 
end 

;=========================================

; Get a radius for each bin.
undef("getRad")
procedure getRad(history:file, shortName:string, X:float, nBins:integer)
  
;  X( ,b) = binned quantity from history (i.e. mmr or radius) 

local b, binName, nBins, xbin

begin
  
  ; Iterate over the bins
  
  do b = 0, nBins-1
      binName = shortName + sprinti("%2.2i", b+1) + "RA"

      ; Get the radius for this bin.
      xbin = history->$binName$

      ; Set array of radii.
      X(:, :, :, :, b) = xbin

      delete(xbin)
      print(binName+" max="+max(X(:, :, :, :, b))+" min="+min(X(:, :, :, :, b)))
;      if (binBame.eq."MIXSUL01RA") then
      if (b.eq.0) then
	a1D      = ndtooned(X(:, :, :, :, b))
	dsizes_a = dimsizes(X)
	indices  = ind_resolve(ind(a1D.eq.max(X(:, :, :, :, b))),dsizes_a)
        print("max ind:"+indices)
	print("time:"+X&time((/indices(0,0)/)))
	print(" lev:"+X&lev((/indices(0,1)/)))
	print(" lat:"+X&lat((/indices(0,2)/)))
	print(" lon:"+X&lon((/indices(0,3)/)))
      end if
  end do
 
end 

;=========================================

; Get a binned quantity for each bin.
undef("getXbin")
procedure getXbin(history:file, shortName:string, X:double, nBins:integer)
  
;  X( ,b) = binned quantity from history (i.e. mmr or radius) 

local b, binName, xbin

begin
  
  print("Getting "+shortName)
  ; Iterate over the bins
  
  do b = 0, nBins-1
      binName = shortName + sprinti("%2.2i", b+1)

      ; Get the quantity for this bin.
      xbin = history->$binName$
;      print("xbin:"+dimsizes(xbin))
;      print("X:"+dimsizes(X))

      ; Set array of binned quantities.
      X(:, :, :, :, b) = xbin

      delete(xbin)
      print(binName+" max="+sprintf("%9.3e",max(X(:, :, :, :, b)))+\
      " min="+sprintf("%9.3e",min(X(:, :, :, :, b))))
  end do
 
end 

;=========================================

; Set a binned quantity for each bin in a history file
undef("setXbin")
procedure setXbin(history:file, shortName:string, X:double, nBins:integer)
  
;  X( ,b) = binned quantity (i.e. mmr or radius) 

local b, binName, xbin

begin
  
  print("Setting "+shortName)
  ; Iterate over the bins
  
  do b = 0, nBins-1
      binName = shortName + sprinti("%2.2i", b+1)
      
      print(binName+" max="+sprintf("%9.3e",max(X(:, :, :, :, b)))+\
           " min="+sprintf("%9.3e",min(X(:, :, :, :, b))))

      ; Get the quantity for this bin.
      xbin = X(:, :, :, :, b)

      ; Set the value in the history file.
      history->$binName$ = xbin

      delete(xbin)
;      print(" bin "+sprinti("%2.2i",b)+" max="+sprintf("%9.3e",max(X(:, :, :, :, b)))+\
;           " min="+sprintf("%9.3e",min(X(:, :, :, :, b))))
  end do
 
end 

;=========================================

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

;=========================================

undef("calcSulfRad")
procedure calcSulfRad(history:file, mixsulnd:float, wtpct:float, sulfDens:float, \
     coreName:string, rmass:float, arho:float, sulfRad:float, coreRad:float, \
     pureRad:float, skipbin:integer)

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
      if (.not.(isfilevar(history,"SULFRAD")).or..not.(isfilevar(history,"CORERAD"))) then
        binName = coreName + sprinti("%2.2i", b+1)

        ; Get the mass mixing ratio.
        mmr = history->$binName$

        ; Compute core mass (g/particle)
        coremass = (mmr * arho) / (mixsulnd (:, :, :, :, b)>1e-30)
        corevol = coremass/RHO_DUST
        coreRad(:, :, :, :, b) = (corevol/cpi)^(1./3.)
        print("bin "+b+", coreRad max="+max(coreRad(:, :, :, :, b))+" min="+min(coreRad(:, :, :, :, b)))

        if (.not.(isfilevar(history,"SULFRAD"))) then
          voltot = corevol + (rmass(b)-coremass) * 100.0/wtpct / sulfDens

          sulfRad(:, :, :, :, b) = (voltot/cpi)^(1./3.)
          print("bin "+b+", sulfRad max="+max(sulfRad(:, :, :, :, b))+" min="+min(sulfRad(:, :, :, :, b)))

          delete(mmr)
        end if
      end if
      
      if (.not.(isfilevar(history,"PURERAD"))) then
        voltot = rmass(b) * 100.0/wtpct / sulfDens

        pureRad(:, :, :, :, b) = (voltot/cpi)^(1./3.)
        print("bin "+b+", pureRad max="+max(pureRad(:, :, :, :, b))+" min="+min(pureRad(:, :, :, :, b)))      
      end if     
      
    end if
    
    if (skipbin.gt.0) then 
      sulfRad(:, :, :, :, skipbin) = sqrt(sulfRad(:, :, :, :, skipbin-1)*sulfRad(:, :, :, :, skipbin+1))
      coreRad(:, :, :, :, skipbin) = sqrt(coreRad(:, :, :, :, skipbin-1)*coreRad(:, :, :, :, skipbin+1))
      pureRad(:, :, :, :, skipbin) = sqrt(pureRad(:, :, :, :, skipbin-1)*pureRad(:, :, :, :, skipbin+1))
    end if
    
  end do
  
end

;=========================================
