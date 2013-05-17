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
  hyai 	= history->hyai
  hybi 	= history->hybi
  P0 	= history->P0
  PS 	= history->PS
  print("  Reading so4 variables ...")
  so4_a1= history->so4_a1        ; kg/kg air
  so4_a2= history->so4_a2        ; kg/kg air
  so4_a3= history->so4_a3        ; kg/kg air
  print("  Summing so4 variables ...")
  so4 = so4_a1 + so4_a2 + so4_a3 ; kg/kg air
  print("  Reading dgnd variables ...")
  dgnd_a1= history->dgnd_a1      ; diameter (m)
  dgnd_a2= history->dgnd_a2      ; diameter (m)
  dgnd_a3= history->dgnd_a3      ; diameter (m)
  print("  Reading num variables ...")
  num_a1= history->num_a1        ; kg/kg air
  num_a2= history->num_a2        ; kg/kg air
  num_a3= history->num_a3        ; kg/kg air
  

  ; Calculate the air pressures at the interfaces.
  print("  Calculating air pressure ...")
  P = pres_hybrid_ccm(PS(:,:,:), P0, hyai(:) ,hybi(:))
  delete(PS)
;  pdel = P(:, 1:nLevs, :, :) - P(:, 0:nLevs-1, :, :)
  
  ; Calculate the atmospheric density (g/cm3) assuming an ideal gas.
  print("  Calculating air density ...")
  arhomass(:,:,:,:)  = P(:,:,:,:) / (RDAIR * T(:,:,:,:))
  arhomass@units = "g/cm3"
  
  so4_nd = so4 * arhomass * AVOGAD / MWNH4HSO4
  so4_nd@units = "molec/cm3"

  
  
end

