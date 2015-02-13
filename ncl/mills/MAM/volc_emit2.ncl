;volc_emit2.ncl
; based on aero2010alt3.ncl
;andrew gettelman nov2014: for volcanic emissions
; aerosol data...
;version history:
;volc_emit2.ncl: add SO4 emissions, and SO4 number concentration additions

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin


fillval=-99999.

;write emissions?
ewrite = "Y"

ipth="/datalocal/haystack/andrew/aero_recent/trop_mozart_aero/emis/"


filedate = systemfunc("date +%Y%m%d")
etxt="Holuhraun35"
;cleanup...daily files and ALT files...

cmd="rm -f "+ipth+"*c"+filedate
system(cmd)
cmd="rm -f "+ipth+"*"+etxt+".nc"
system(cmd)

;sector names
sectxt=(/"awb","tra","wst","dom","ind","ene"/)
ns=dimsizes(sectxt)
vars="emis_"+sectxt

;------------------------------------------------------------
;PARAMETERS
;------------------------------------------------------------

;Holuhran version 2 (corrected)
;if(etxt.eq."Holuhraun3") then 
  vlat=63.7
  vlon=360-16.8

;	Have altitude intervals: ...1.398, 1.884, 2.579, 3.471...
; altitude = 0.063, 0.202, 0.366, 0.554, 0.767, 1.003, 1.262, 1.641, 2.232, 
;    3.025, 4.009, 5.157, 6.356 ;
;altitude_int = 0, 0.126, 0.278, 0.454, 0.654, 0.879, 1.127, 1.398, 1.884, 
;    2.579, 3.471, 4.548, 5.765, 6.946 ;
;		(so goes into 3 grid boxes)  
;	(m in vertical)  486m,695m,892m  = 4.86e4cm, 6.95e4cm, 8.92e4cm

;be careful of levels. Watch intervals

if (etxt.eq."HoluhraunHi") then                 ;Hi 3-4.5km
  altcm=(/892,1077./)*100.
  valt=(/2.99,4.5/)  ;intefaces 2.579-4.548
else if (etxt.eq."HoluhraunLo") then            ;LO: 0.5-2km
  altcm=(/200.,225,248,271,486. /)*100. 
  valt=(/0.5,2.0/)  ;interfaces 0.454-1.884
else                                     ;BASE 1.5-3km
  altcm=(/486,695,892./)*100.
  valt=(/1.5,3.2/)  ;interfaces 1.398-3.471
end if
end if 

print(altcm)
print(valt)



;so2 emissions in molecules/s
  so2ktd=35.
  so2kgd=so2ktd*1000.*1000.

  nakg=6.022e26
  secday=3600.*24.
  molmass= 64.

  so2molecs=so2kgd * nakg/secday/64.

print(so2molecs)

;At 63.7 N, one grid box is:
  dx=2.5
  dy=1.9
  re=6371.
  pi=3.14159

  dykm= dy * re * 2.*pi/360. 
  dxkm=	dx * re * 2 * pi/360. * cos(2*pi*vlat/360.) 

  ;print(dykm)
  ;print(dxkm)

  km2cm=1000.*100.

  area=(dykm*km2cm)*(dxkm*km2cm)
  ;print(area)

  vol=area*sum(altcm)

; now convert...

  so2val=so2molecs/vol
  print(so2val)

;  so2val=5.233e9 INCORRECT Value (version 1 factor of 64 too high!)
;  so2val=8.178e7 ;corrected for molecular mass...molecules/cm3/s
;  now calculate online as above from kT SO2/day

;end if
;EMISSIONS FRACTIONS INTO SO2 and parameters
; SO4 aiken = a2
; SO4 accum = a1

; to convert g cm-3 to kg m-3  * 1.e6 cm3 m-3 / 1.e3 g kg-1
densconv = 1.e6 / 1.e3  

fracso2=0.975

; emissions fraction is multiplied by 1.80 to get ammonium-bisulfate mass
fmult=1.80

pi=3.14159

;a1
fso4a1=0.025 * fmult
rhoso4a1= 1.77 * densconv ;density (g/cm-3, convert to kg)
dso4a1= 0.134  * 1.e-6 ; diameter (microns, convert to m)

;mass/particle (kg)
a1pmass = rhoso4a1 * 4./3. * pi * (dso4a1/2.)^3

;conversion: want molecules/particle  * 6.022e26 

;kg/particle *  (molecules/kg mol) / kg/kg mol / (molecules/kg mol)
;check where the molecular mass goes...? If at all?
;heavier masses make for fewer molecules per particle?  

molp1= a1pmass / 96.

;a2  NO A2 ELEVATED EMISSIONS
;fso4a2=0.0125 * fmult
;rhoso4a2= 1.77 * densconv ;density (g/cm-3, convert to kg)
;dso4a2= 0.0504  * 1.e-6 ; diameter (microns, convert to m)

;;mass/particle (kg)
;a2pmass = rhoso4a2 * 4./3. * pi * (dso4a2/2.)^3

;conversion: want molecules/particle  * 6.022e26 

;kg/particle *  (molecules/kg mol) / kg/kg mol / (molecules/kg mol)
;check where the molecular mass goes...? If at all?
;heavier masses make for fewer molecules per particle?  

;molp2= a2pmass / 96.

;------------------------------------------------------------
;READ data
;------------------------------------------------------------
;originally from:  /fs/cgd/csm/inputdata/atm/cam/chem/trop_mozart_aero/emis/

so2_elev_2000="ar5_mam3_so2_elev_2000_c090726.nc"

ie2000=addfile(ipth+so2_elev_2000,"r")

;read base emissions
ev= ie2000->contvolc

ev2=ev
;printVarSummary(ev)
printVarSummary(ev2)

;print(ev(0,{1.5:3.2},{vlat},{vlon}))

;mass of SO2 = molecules/cm3/s

ev2(:,{valt(0):valt(1)},{vlat},{vlon})=ev2(:,{valt(0):valt(1)},{vlat},{vlon})+so2val*fracso2

print(ev2(0,{valt(0):valt(1)},{vlat},{vlon}))

;------------------------------------------------------------
;SO4 Mass
;------------------------------------------------------------

;------------------------------------------------------------
;a1
so4a1_elev_2000="ar5_mam3_so4_a1_elev_2000_c090726.nc"

i41e2000=addfile(ipth+so4a1_elev_2000,"r")

;read base emissions
e41v= i41e2000->contvolc

e41v2=e41v

;printVarSummary(e41v2)

;mass of SO4 = molecules/cm3/s

e41v2(:,{valt(0):valt(1)},{vlat},{vlon})=e41v2(:,{valt(0):valt(1)},{vlat},{vlon})+so2val*fso4a1

;------------------------------------------------------------
;a2; surface only (neglect for now?)

;------------------------------------------------------------
;SO4 Number
;------------------------------------------------------------

;------------------------------------------------------------
;a1
na1_elev_2000="ar5_mam3_num_a1_elev_2000_c090726.nc"

in1e2000=addfile(ipth+na1_elev_2000,"r")

;read base emissions
n1v= in1e2000->SO4_contvolc

n1v2=n1v

;printVarSummary(e41v2)

;a1 number = particles/cm3/s = molecules/cm3/s  / molecules/parrticle

n1v2(:,{valt(0):valt(1)},{vlat},{vlon})=n1v2(:,{valt(0):valt(1)},{vlat},{vlon})+so2val*fso4a1 /molp1

;------------------------------------------------------------
;test global mean....
;------------------------------------------------------------
;;weight
;kwgty=NormCosWgtGlobe(skli2000&lat)/2.
;kwgty!0="lat"
;kwgty&lat=skli2000&lat

;rwgty=NormCosWgtGlobe(f2000&lat)/2.
;rwgty!0="lat"
;rwgty&lat=f2000&lat

;;printVarSummary(wgty)

;;global means
;gmk2000 = wgt_areaave_Wrap(skli2000,kwgty,1.0,0)
;gms2000 = wgt_areaave_Wrap(f2000,rwgty,1.0,0)

;print("----------------")
;print("Weighted annual avg")
;print("Original Grid ="+gmk2000)
;print("Regridded="+gms2000)
;print("Percent Difference="+(gms2000-gmk2000)/gmk2000*100.)
;print("----------------")



;------------------------------------------------------------
;now, read volcanic emissions, process and write back into revised file....
;------------------------------------------------------------

;sector maps...

iawb=ind(sectxt.eq."awb")
itra=ind(sectxt.eq."tra")
iind=ind(sectxt.eq."ind")
iwst=ind(sectxt.eq."wst")
idom=ind(sectxt.eq."dom")
iene=ind(sectxt.eq."ene")
;print(idom)
;exit

if (ewrite.eq."Y") then
;------------------------------------------------------------
;s02 elev

;write to revised file...
  outSO2="ar5_mam3_so2_elev_"+etxt+"_c"+filedate+".nc"
  system("cp "+ipth+so2_elev_2000+" "+ipth+outSO2)

  so2s=addfile(ipth+outSO2,"w")

  so2s->contvolc =ev2

  delete(so2s)


;------------------------------------------------------------
;s04 a1 elev

;write to revised file...
  outSO41="ar5_mam3_so4_a1_elev_"+etxt+"_c"+filedate+".nc"
  system("cp "+ipth+so4a1_elev_2000+" "+ipth+outSO41)

  so4s1=addfile(ipth+outSO41,"w")

  so4s1->contvolc =e41v2

  delete(so4s1)

;------------------------------------------------------------
;num a1 elev

;write to revised file...
  outn1="ar5_mam3_num_a1_elev_"+etxt+"_c"+filedate+".nc"
  system("cp "+ipth+na1_elev_2000+" "+ipth+outn1)

  so4n1=addfile(ipth+outn1,"w")

  so4n1->SO4_contvolc =n1v2

  delete(so4n1)


end if ;end write conditional





end
