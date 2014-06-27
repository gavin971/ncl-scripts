load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

MWH2SO4 = 98.0784                 ; molecular weight of H2SO4 species (g/mol)
MWS     = 32.066                  ; molecular weight of atomic sulfur (g/mol)

caseid=getenv("caseid")

Rearth  = 6.37122e6               ; Earth's radius (m)
PI      = 4.0 * atan(1.0)         ; pi

SAearth=4.e0*PI*Rearth*Rearth ; Earth surface area in m2 (CAM)
;print("SAearth="+SAearth+" m2")
cfac = SAearth/1.e6 ; (kgS/m2 -> TgS)
;print(cfac)


if (ismissing(caseid)) then
  print("missing environment variable caseid, quitting")
else
  rootpath=getenv("rootpath")
  if (ismissing(rootpath)) then
    rootpath="/glade/scratch/mmills/archive/"
  end if
  gmdir=rootpath+"/"+caseid+"/atm/proc/h0gm/"
  filename=caseid+".cam.h0.globavg.nc"
  a=addfile(gmdir+filename,"r")
  TMso4_strat   = a->TMso4_strat
  printVarSummary(TMso4_strat)
  time=TMso4_strat&time
  year1=(time(0)-31)/365
  nMonths=dimsizes(TMso4_strat)
  nYears=nMonths/12
  TMso4_strat   = (/ TMso4_strat * cfac /)
  if (isfilevar(a,"TMso4_trop")) then
    TMso4_trop    = a->TMso4_trop
    TMSO2_strat   = a->TMSO2_strat
    TMSO2_trop    = a->TMSO2_trop
    TMOCS_strat   = a->TMOCS_strat
    TMOCS_trop    = a->TMOCS_trop
    TMDMS_strat   = a->TMDMS_strat
    TMDMS_trop    = a->TMDMS_trop
    TMH2SO4_strat = a->TMH2SO4_strat
    TMH2SO4_trop  = a->TMH2SO4_trop

    TMso4_trop    = (/ TMso4_trop * cfac /)
    TMSO2_strat   = (/ TMSO2_strat * cfac /)
    TMSO2_trop    = (/ TMSO2_trop * cfac /)
    TMOCS_strat   = (/ TMOCS_strat * cfac /)
    TMOCS_trop    = (/ TMOCS_trop * cfac /)
    TMDMS_strat   = (/ TMDMS_strat * cfac /)
    TMDMS_trop    = (/ TMDMS_trop * cfac /)
    TMH2SO4_strat = (/ TMH2SO4_strat * cfac /)
    TMH2SO4_trop  = (/ TMH2SO4_trop * cfac /)

    print((/caseid/))
    print("time(0)="+time(0)+" "+TMso4_strat&time@units)
    print("Stratospheric burdens (Tg S):")
;    print("OCS   DMS      SO2  H2SO4 sulfate")
;    print(sprintf("%5.1f",TMOCS_strat)+" "+sprintf("%5.1e",TMDMS_strat)+" "+sprintf("%5.2f",TMSO2_strat)+" "+\
;  	  sprintf("%5.3f",TMH2SO4_strat)+" "+sprintf("%5.1f",TMso4_strat))
    print("Year OCS   DMS      SO2  H2SO4 sulfate")
    do i = 0,nYears-1
      t1=12*i
      t2=t1+11
      year = i+1
      print(" "+year +"   " + sprintf("%5.1f",avg(TMOCS_strat(t1:t2)))+" "+\
            sprintf("%5.1e",avg(TMDMS_strat(t1:t2)))+" "+sprintf("%5.2f",avg(TMSO2_strat(t1:t2)))+" "+\
  	    sprintf("%5.3f",avg(TMH2SO4_strat(t1:t2)))+" "+sprintf("%5.1f",avg(TMso4_strat(t1:t2))))
    end do
    
    print("Tropospheric burdens (Tg S):")
;    print("OCS    DMS    SO2  H2SO4 sulfate")
;    print(sprintf("%5.1f",TMOCS_trop)+" "+sprintf("%5.1f",TMDMS_trop)+" "+sprintf("%5.2f",TMSO2_trop)+" "+\
;  	  sprintf("%5.3f",TMH2SO4_trop)+" "+sprintf("%5.1f",TMso4_trop))
    print("Year OCS   DMS      SO2  H2SO4 sulfate")
    do i = 0,3
      t1=12*i
      t2=t1+11
      year = i+1
      print(" "+year +"   " + sprintf("%5.1f",avg(TMOCS_trop(t1:t2)))+" "+\
            sprintf("%5.1f",avg(TMDMS_trop(t1:t2)))+" "+sprintf("%5.2f",avg(TMSO2_trop(t1:t2)))+" "+\
  	    sprintf("%5.3f",avg(TMH2SO4_trop(t1:t2)))+" "+sprintf("%5.1f",avg(TMso4_trop(t1:t2))))
    end do
  else
    TMso4_strat=TMso4_strat*MWS/MWH2SO4
    print((/caseid/))
    print("time(0)="+time(0)+" "+TMso4_strat&time@units)
    print( "Stratospheric sulfate burdens (Tg S):" )
;    print( sprintf("%5.1f",TMso4_strat) )

    print("Year sulfate")
    do i = 0,nYears-1
      t1=12*i
      t2=t1+11
      year = i+1
      print( " "+year +"   " + sprintf("%5.1f",avg(TMso4_strat(t1:t2))) )   
    end do
  end if
end if
