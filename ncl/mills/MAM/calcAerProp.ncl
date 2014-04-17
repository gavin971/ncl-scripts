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
MWNH4HSO4 = 115.108999            ; molecular weight for MAM so4 species (ammonium sulfate)
MWH2SO4 = 98.078479               ; molecular weight of H2SO4 species

procedure calcAerProp(rootpath:string,filename:string,sigma[3]:float)

begin

  print("reading "+filename)
  history=addfile(rootpath+filename,"w")
  
  ; Compute the air mass and density.
  ;
  ; NOTE: Convert arho and amass are in cgs units.
  print("")
  print("  Reading pressure variables ...")
  hyam 	= history->hyam
  hybm 	= history->hybm
  hyai 	= history->hyai
  hybi 	= history->hybi
  P0 	= history->P0
  PS 	= history->PS
  
  if (isfilevar(history,"so4")) then
    print("  Reading so4 ...")
    so4 = history->so4
    printVarSummary(so4)
  else
    print("  Reading so4 variables ...")
    so4_a1= history->so4_a1        ; kg/kg air
    so4_a2= history->so4_a2        ; kg/kg air
    so4_a3= history->so4_a3        ; kg/kg air
    print("  Summing so4 variables ...")
    so4 = so4_a1 ; get metadata
    so4 = so4_a1 + so4_a2 + so4_a3 ; kg/kg air
    so4@long_name="so4 concentration"
    print("  Saving so4 ...")
    printVarSummary(so4)
    history->so4 = so4
  end if
  
  if (isfilevar(history,"TMso4")) then
    print("  Reading TMso4 ...")
    TMso4=history->TMso4
  else
    print("  Reading TMso4 variables ...")
    TMso4_a1= history->TMso4_a1        ; kg/kg air
    TMso4_a2= history->TMso4_a2        ; kg/kg air
    TMso4_a3= history->TMso4_a3        ; kg/kg air
    print("  Summing TMso4 variables ...")
    TMso4 = TMso4_a1 ; get metadata
    TMso4 = TMso4_a1 + TMso4_a2 + TMso4_a3 ; kg/kg air
    TMso4@long_name="so4 column burden"
    print("  Saving TMso4 ...")
    printVarSummary(TMso4)
    history->TMso4 = TMso4
  end if
  
;  print("  Reading dgnd variables ...")
;  dgnd_a1= history->dgnd_a1      ; diameter (m)
;  dgnd_a2= history->dgnd_a2      ; diameter (m)
;  dgnd_a3= history->dgnd_a3      ; diameter (m)
;  print("  Reading num variables ...")
;  num_a1= history->num_a1        ; kg/kg air
;  num_a2= history->num_a2        ; kg/kg air
;  num_a3= history->num_a3        ; kg/kg air
  

  ; Calculate the air pressures at the interfaces.
  if (isfilevar(history,"P")) then
    print("  Reading P ...")
    P = history->P
    history->P@units="Pa"
    P@units="Pa"
  else
    print("  Calculating air pressure ...")
    P = so4 ; get metadata
    P = pres_hybrid_ccm(PS(:,:,:), P0, hyam(:) ,hybm(:))
    P@long_name="Pressure"
    P@units="Pa"   
    print("  Saving P ...")
    printVarSummary(P)
    history->P = P
  end if
  
  if (isfilevar(history,"TMso4_strat")) then
    print("  Found TMso4_strat ...")
;    TMso4_strat = history->TMso4_strat
;    TMso4_strat = TMso4_strat*MWH2SO4/MWNH4HSO4 ; scale molecular weight
;    TMso4_strat@method="scaled by 98.078479/115.108999 to adjust mol wt NH4HSO4 to H2SO4"
;    print("  Saving TMso4_strat ...")
;    printVarSummary(TMso4_strat)
;    history->TMso4_strat = TMso4_strat
  else
    TMso4_strat=TMso4 ; get metadata and array structure
    print("  Reading TROP_P ...")
    tropP= history->TROP_P        ; tropopause pressure
    print("calculating pdel ...")
    Pi = pres_hybrid_ccm(PS(:,:,:), P0, hyai(:) ,hybi(:))   
    nLevs=dimsizes(so4&lev)
    pdel = Pi(:, 1:nLevs, :, :) - Pi(:, 0:nLevs-1, :, :) ; (time, lev, lat, lon)    
    cellMass = so4 ; get metadata
    cellMass = (pdel / G) * so4 ; k/m2 per level
    nLat=dimsizes(tropP&lat)
    nLon=dimsizes(tropP&lon)
    do j = 0,nLat-1
      do m = 0,nLon-1
        tp = tropP(0,j,m)
        TMso4_strat(0,j,m) = dim_sum(cellMass(time|0, lat|j, lon|m, {lev|:tp}))
      end do
    end do    
    delete(cellMass)
    TMso4_strat=TMso4_strat*MWH2SO4/MWNH4HSO4

    TMso4_strat@long_name="so4 stratospheric column burden"
    TMso4_strat@method="scaled by 98.078479/115.108999 to adjust mol wt NH4HSO4 to H2SO4"
    print("  Saving TMso4_strat ...")
    printVarSummary(TMso4_strat)
    history->TMso4_strat = TMso4_strat
  end if
  
  ; Calculate the atmospheric density (g/cm3) assuming an ideal gas.
  if (isfilevar(history,"arhomass")) then
    print("  Reading arhomass ...")
    arhomass = history->arhomass
    printVarSummary(arhomass)    
  else
    print("  Reading T ...")
    T = history->T
    printVarSummary(T)
    arhomass = T ; get metadata
    print("  Calculating air density ...")
    arhomass(:,:,:,:)  = P(:,:,:,:) / (RDAIR * T(:,:,:,:)) / 1.e3 ; kg/m3 -> g/cm3
    arhomass@long_name="air density"
    arhomass@units = "g/cm3"
    print("  Saving arhomass ...")
    printVarSummary(arhomass)    
    history->arhomass=arhomass
  end if
  
  if (isfilevar(history,"so4_nd")) then
    print("  Found so4_nd ...")
;    print("  Reading so4_nd ...")
;    so4_nd = history->so4_nd
;    printVarSummary(so4_nd)    
  else
    print("  Calculating sulfate molec/cm3 ...")
    so4_nd = so4 ; get metadata
    so4_nd = so4(:,:,:,:) * arhomass(:,:,:,:) * AVOGAD / 1.e3 / MWNH4HSO4
    so4_nd@long_name="so4 number density"
    so4_nd@units = "molec/cm3"
    print("  Saving so4_nd ...")
    printVarSummary(so4_nd)    
    history->so4_nd=so4_nd
  end if

;  if (isfilevar(history,"radeffwet")) then
;    print("  Found radeffwet ...")
;  else
      
    print("  Reading dgnumwet1 ...")
    diam1=history->dgnumwet1
    print("  Reading dgnumwet2 ...")
    diam2=history->dgnumwet2
    print("  Reading dgnumwet3 ...")
    diam3=history->dgnumwet3
    
    exists=isvar((/"num_a1","num_a2","num_a3"/))
    if(.not.exists(0)) then
      print("  Reading num_a1 ...")
      num_a1= history->num_a1        ; kg/kg air
    end if
    if(.not.exists(1)) then
      print("  Reading num_a2 ...")
      num_a2= history->num_a2        ; kg/kg air
    end if
    if(.not.exists(2)) then
      print("  Reading num_a3 ...")
      num_a3= history->num_a3        ; kg/kg air
    end if
        
    print("  Calculating surfwet ...")
    surf_a1=diam1 ; get metadata
    surf_a1@units="um2/cm3"
    surf_a2=surf_a1 ; get metadata
    surf_a3=surf_a1 ; get metadata
    surf_a1=num_a1*diam1^2.*exp(2.*log(sigma(0))^2.)
    surf_a2=num_a2*diam2^2.*exp(2.*log(sigma(1))^2.)
    surf_a3=num_a3*diam3^2.*exp(2.*log(sigma(2))^2.)
    surf_a1=surf_a1*1e12*1e-3*arhomass*PI ; convert m2/kg air -> um2/cm3 air
    surf_a2=surf_a2*1e12*1e-3*arhomass*PI ; convert m2/kg air -> um2/cm3 air
    surf_a3=surf_a3*1e12*1e-3*arhomass*PI ; convert m2/kg air -> um2/cm3 air
    surf = surf_a1  ; get metadata
    surf = surf_a1  + surf_a2 + surf_a3
    surf@_FillValue=9.96921e+36
    surf=where(surf.gt.0,surf,surf@_FillValue)
    surf@long_name="Aerosol wet surface area"
    surf_a1@_FillValue=9.96921e+36
    surf_a1=where(surf_a1.gt.0,surf_a1,surf_a1@_FillValue)
    surf_a1@long_name="Accumulation mode wet surface area"
    surf_a2@_FillValue=9.96921e+36
    surf_a2=where(surf_a2.gt.0,surf_a2,surf_a2@_FillValue)
    surf_a2@long_name="Nucleation mode wet surface area"
    surf_a3@_FillValue=9.96921e+36
    surf_a3=where(surf_a3.gt.0,surf_a3,surf_a3@_FillValue)
    surf_a3@long_name="Coarse mode wet surface area"
   
    print("  Saving surfwet ...")
    printVarSummary(surf)    
    history->surfwet=surf
    print("  Saving surfwet_a1 ...")
    printVarSummary(surf_a1)    
    history->surfwet_a1=surf_a1
    print("  Saving surfwet_a2 ...")
    printVarSummary(surf_a2)    
    history->surfwet_a2=surf_a2
    print("  Saving surfwet_a3 ...")
    printVarSummary(surf_a3)    
    history->surfwet_a3=surf_a3
    
    print("  Calculating volwet ...")
    vol_a1=diam1 ; get metadata
    vol_a1@units="um3/cm3"
    vol_a2=vol_a1 ; get metadata
    vol_a3=vol_a1 ; get metadata
    vol_a1=num_a1*diam1^3.*exp(3.*log(sigma(0))^2.)
    vol_a2=num_a2*diam2^3.*exp(3.*log(sigma(1))^2.)
    vol_a3=num_a3*diam3^3.*exp(3.*log(sigma(2))^2.)
    vol_a1=vol_a1*1e18*1e-3*arhomass*PI/6. ; convert m3/kg air -> um3/cm3 air
    vol_a2=vol_a2*1e18*1e-3*arhomass*PI/6. ; convert m3/kg air -> um3/cm3 air
    vol_a3=vol_a3*1e18*1e-3*arhomass*PI/6. ; convert m3/kg air -> um3/cm3 air
    vol = vol_a1  ; get metadata
    vol = vol_a1  + vol_a2 + vol_a3
    vol@_FillValue=9.96921e+36
    vol=where(vol.gt.0,vol,vol@_FillValue)
    vol@long_name="Aerosol wet volume"
    vol_a1@_FillValue=9.96921e+36
    vol_a1=where(vol_a1.gt.0,vol_a1,vol_a1@_FillValue)
    vol_a1@long_name="Accumulation mode wet volume"
    vol_a2@_FillValue=9.96921e+36
    vol_a2=where(vol_a2.gt.0,vol_a2,vol_a2@_FillValue)
    vol_a2@long_name="Nucleation mode wet volume"
    vol_a3@_FillValue=9.96921e+36
    vol_a3=where(vol_a3.gt.0,vol_a3,vol_a3@_FillValue)
    vol_a3@long_name="Coarse mode wet volume"
   
    print("  Saving volwet ...")
    printVarSummary(vol)    
    history->volwet=vol
    print("  Saving volwet_a1 ...")
    printVarSummary(vol_a1)    
    history->volwet_a1=vol_a1
    print("  Saving volwet_a2 ...")
    printVarSummary(vol_a2)    
    history->volwet_a2=vol_a2
    print("  Saving volwet_a3 ...")
    printVarSummary(vol_a3)    
    history->volwet_a3=vol_a3

    print("  Calculating radeffwet ...")
    radeffwet=diam1 ; get metadata
    radeffwet=3.0*vol/surf
    radeffwet@long_name="Aerosol wet effective radius"
    radeffwet@units="um"
   
    print("  Saving radeffwet ...")
    printVarSummary(radeffwet)    
    history->radeffwet=radeffwet

;  end if

;  if (isfilevar(history,"radeffdry")) then
;    print("  Found radeffdry ...")
;  else
      
    delete(diam1)
    delete(diam2)
    delete(diam3)
    print("  Reading dgnd_a01 ...")
    diam1=history->dgnd_a01
    print("  Reading dgnd_a02 ...")
    diam2=history->dgnd_a02
    print("  Reading dgnd_a03 ...")
    diam3=history->dgnd_a03
    
    exists=isvar((/"num_a1","num_a2","num_a3"/))
    if(.not.exists(0)) then
      print("  Reading num_a1 ...")
      num_a1= history->num_a1        ; kg/kg air
    end if
    if(.not.exists(1)) then
      print("  Reading num_a2 ...")
      num_a2= history->num_a2        ; kg/kg air
    end if
    if(.not.exists(2)) then
      print("  Reading num_a3 ...")
      num_a3= history->num_a3        ; kg/kg air
    end if
        
    delete(surf)
    print("  Calculating surfdry ...")
    surf_a1=diam1 ; get metadata
    surf_a1@units="um2/cm3"
    surf_a2=surf_a1 ; get metadata
    surf_a3=surf_a1 ; get metadata
    surf_a1=num_a1*diam1^2.*exp(2.*log(sigma(0))^2.)
    surf_a2=num_a2*diam2^2.*exp(2.*log(sigma(1))^2.)
    surf_a3=num_a3*diam3^2.*exp(2.*log(sigma(2))^2.)
    surf_a1=surf_a1*1e12*1e-3*arhomass*PI ; convert m2/kg air -> um2/cm3 air
    surf_a2=surf_a2*1e12*1e-3*arhomass*PI ; convert m2/kg air -> um2/cm3 air
    surf_a3=surf_a3*1e12*1e-3*arhomass*PI ; convert m2/kg air -> um2/cm3 air
    surf = surf_a1  ; get metadata
    surf = surf_a1  + surf_a2 + surf_a3
    surf@_FillValue=9.96921e+36
    surf=where(surf.gt.0,surf,surf@_FillValue)
    surf@long_name="Aerosol dry surface area"
    surf_a1@_FillValue=9.96921e+36
    surf_a1=where(surf_a1.gt.0,surf_a1,surf_a1@_FillValue)
    surf_a1@long_name="Accumulation mode dry surface area"
    surf_a2@_FillValue=9.96921e+36
    surf_a2=where(surf_a2.gt.0,surf_a2,surf_a2@_FillValue)
    surf_a2@long_name="Nucleation mode dry surface area"
    surf_a3@_FillValue=9.96921e+36
    surf_a3=where(surf_a3.gt.0,surf_a3,surf_a3@_FillValue)
    surf_a3@long_name="Coarse mode dry surface area"
   
    print("  Saving surfdry ...")
    printVarSummary(surf)    
    history->surfdry=surf
    print("  Saving surfdry_a1 ...")
    printVarSummary(surf_a1)    
    history->surfdry_a1=surf_a1
    print("  Saving surfdry_a2 ...")
    printVarSummary(surf_a2)    
    history->surfdry_a2=surf_a2
    print("  Saving surfdry_a3 ...")
    printVarSummary(surf_a3)    
    history->surfdry_a3=surf_a3
    
    delete(vol)
    print("  Calculating voldry ...")
    vol_a1=diam1 ; get metadata
    vol_a1@units="um3/cm3"
    vol_a2=vol_a1 ; get metadata
    vol_a3=vol_a1 ; get metadata
    vol_a1=num_a1*diam1^3.*exp(3.*log(sigma(0))^2.)
    vol_a2=num_a2*diam2^3.*exp(3.*log(sigma(1))^2.)
    vol_a3=num_a3*diam3^3.*exp(3.*log(sigma(2))^2.)
    vol_a1=vol_a1*1e18*1e-3*arhomass*PI/6. ; convert m3/kg air -> um3/cm3 air
    vol_a2=vol_a2*1e18*1e-3*arhomass*PI/6. ; convert m3/kg air -> um3/cm3 air
    vol_a3=vol_a3*1e18*1e-3*arhomass*PI/6. ; convert m3/kg air -> um3/cm3 air
    vol = vol_a1  ; get metadata
    vol = vol_a1  + vol_a2 + vol_a3
    vol@_FillValue=9.96921e+36
    vol=where(vol.gt.0,vol,vol@_FillValue)
    vol@long_name="Aerosol dry volume"
    vol_a1@_FillValue=9.96921e+36
    vol_a1=where(vol_a1.gt.0,vol_a1,vol_a1@_FillValue)
    vol_a1@long_name="Accumulation mode dry volume"
    vol_a2@_FillValue=9.96921e+36
    vol_a2=where(vol_a2.gt.0,vol_a2,vol_a2@_FillValue)
    vol_a2@long_name="Nucleation mode dry volume"
    vol_a3@_FillValue=9.96921e+36
    vol_a3=where(vol_a3.gt.0,vol_a3,vol_a3@_FillValue)
    vol_a3@long_name="Coarse mode dry volume"
   
    print("  Saving voldry ...")
    printVarSummary(vol)    
    history->voldry=vol
    print("  Saving voldry_a1 ...")
    printVarSummary(vol_a1)    
    history->voldry_a1=vol_a1
    print("  Saving voldry_a2 ...")
    printVarSummary(vol_a2)    
    history->voldry_a2=vol_a2
    print("  Saving voldry_a3 ...")
    printVarSummary(vol_a3)    
    history->voldry_a3=vol_a3

    print("  Calculating radeffdry ...")
    radeffdry=diam1 ; get metadata
    radeffdry=3.0*vol/surf
    radeffdry@long_name="Aerosol dry effective radius"
    radeffdry@units="um"
    
    print("  Saving radeffdry ...")
    printVarSummary(radeffdry)    
    history->radeffdry=radeffdry

;  end if
  
end

;caseid="b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001"
;caseid="b.e11.B55TRW5CN.f19_g16.Pin20-30.dgnh1um.noqbo.001"
;caseid="b.e12.BRCP45W5CN.f19_g16.control.002"
;caseid="b.e12.BRCP45W5CN.f19_g16.geo02Tgpy2015TR.002"
;caseid="b.e12.BRCP45W5CN.f19_g16.hifreq.001"
;caseid="b.e12.BRCP45W5CN.f19_g16.geo10Tgpy.001"
;caseid="f.e12.FSDW5.f19_f19.Pin15Tg.19-27km.0N-20N.120E.a3dgnl44.001"
;caseid="f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.a3dgnl44.001"
caseid="f.e12.FSDW5.f19_f19.geo00.02"
;caseid=getenv("caseid")
print(caseid)
casedir=getenv("casedir")
print(casedir)
rootpath="/glade/scratch/mmills/archive/"+caseid+"/atm/proc/h0mam/"
sigma=(/1.6, 1.6, 1.2/)

copypattern="/glade/scratch/mmills/archive/"+caseid+"/atm/hist/*h0*"
if (casedir.eq."run") then
  copypattern="/glade/scratch/mmills/"+caseid+"/run/"+caseid+".cam.h0*"
end if

print("===========================================================")
print("mkdir -p "+rootpath)
mkdir=systemfunc("mkdir -p "+rootpath)
print("===========================================================")
print("cp -n -v "+copypattern+" "+rootpath)
cp=systemfunc("cp -n -v "+copypattern+" "+rootpath)
print(cp)
print("===========================================================")
historyNames=systemfunc("ls "+rootpath)
print(historyNames)
nfiles=dimsizes(historyNames)
print("===========================================================")

do i = 0,nfiles-1
  i1=i+1
  print("===========================================================")
  print("File "+i1+" of "+nfiles+":"+historyNames(i))
  calcAerProp(rootpath,historyNames(i),sigma)
  print("===========================================================")
end do
