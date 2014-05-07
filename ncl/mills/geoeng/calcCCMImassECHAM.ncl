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
REARTH  = 6.37122e6               ; Earth's radius (m)
DG2RAD  = PI / 180.0              ; rad deg^-1

MWS = 32.065 ; molecular weight of S (kg/kmol)

load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;rootpath="/Users/mmills/Dropbox/Work/GeoMIP/"
;filename="background_2001_ccmval_ymon_zm_pl.nc"
rootpath="/glade/p/work/mmills/inputdata/sulf/"
filename="svc002_2001_ccmval_sad_ymon_zm_pl.nc"

a=addfile(rootpath+filename,"w")

lev=a->lev ; pressure (Pa), 39 levels
temp=a->st ; temperature (K), (time=12, lev=39, lat=64, lon=1)
so4=a->SO4_AS
so4_cs=a->SO4_CS
so4=so4+so4_cs ; kg S/kg air

;vol=a->vol_as ; um3/cm2
;vol=vol+a->vol_cs
vol=a->VAD ; m3/m3 = cm3/cm3
vol = where(vol .le. 0.0, vol@_FillValue, vol*1e6) ; cm3/cm3 -> um3/cm3
vol@long_name="aerosol volume"
vol@units="um3/cm3"
vol = where(vol .le. 0.0, vol@_FillValue, vol)
print("saving volume_density")
a->volume_density=vol

;sad=a->aer_surf_as ; um2/cm2
;sad=sad+a->aer_surf_cs
sad=a->SAD ; m2/m3 = cm3/cm3
sad = where(sad .le. 0.0, sad@_FillValue, sad*1e4) ; cm2/cm3 -> um2/cm3
sad@units = "um2/cm3"
print("saving sad")
a->sad=sad

reff = vol ; get metadata
reff = 3*vol/sad ; "um"
reff@long_name="effective radius"
reff@units = "microns"
print("saving reff")
a->reff=reff

; assume mode width 1.25, rmode = reff * exp(-2.5 * (ln(sigma))**2)
rmode =  reff ; get metadata
rmode =  reff* 0.8829536147077603 
rmode@long_name="mode radius"
rmode@method="reff * 0.8829536147077603, assuming mode width 1.25, rmode = reff * exp(-2.5 * (ln(sigma))**2)"
print("saving rmode")
a->rmode=rmode

nlev = dimsizes(lev)

; Calculate the density assuming an ideal gas.
arho = temp ; get metadata
do k = 0,nlev-1
  arho(:,k,:,:)  = dble2flt(lev(k)) / (RDAIR * temp(:,k,:,:)) ; Pa*kg/J = (kg2 m-1 s-2)/(kg m2 s-2) = kg m-3
end do
arho@long_name="air density"
arho@units="kg m-3"

H2SO4_mass = temp ; get metadata
H2SO4_mass = arho*so4 ; kg S m-3
H2SO4_mass = H2SO4_mass*AVOGAD/MWS ; molec m-3
H2SO4_mass = H2SO4_mass*1.e-6 ; molec cm-3
H2SO4_mass@long_name="sulfate aerosol mass"
H2SO4_mass@units="molecules/cm3air"
print("saving H2SO4_mass")
a->H2SO4_mass=H2SO4_mass
