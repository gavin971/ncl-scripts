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
MWNH4HSO4 = 115.108999 

procedure calcAerProp(rootpath:string,filename:string)

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
    delete(PS)
    P@long_name="Pressure"
    P@units="Pa"   
    print("  Saving P ...")
    printVarSummary(P)
    history->P = P
  end if
;  pdel = P(:, 1:nLevs, :, :) - P(:, 0:nLevs-1, :, :)
  
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
    print("  Reading so4_nd ...")
    so4_nd = history->so4_nd
    printVarSummary(so4_nd)    
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

  
end

;caseid="b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001"
;caseid="b.e11.B55TRW5CN.f19_g16.Pin20-30.dgnh1um.noqbo.001"
;caseid="b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnh1um.qboon.001"
caseid=getenv("caseid")
print(caseid)
casedir=getenv("casedir")
print(casedir)
rootpath="/glade/scratch/mmills/archive/"+caseid+"/atm/proc/h0mam/"

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
  calcAerProp(rootpath,historyNames(i))
  print("===========================================================")
end do

