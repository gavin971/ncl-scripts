load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

MWH2SO4 = 98.0784                 ; molecular weight of H2SO4 species (g/mol)
MWS     = 32.066                  ; molecular weight of atomic sulfur (g/mol)

; Global H2SO4 burden (Tg) data from Baran and Foot, JGR, 1994, Figure 3:
BaranMonth = (/ 199105,199106,199107,199108,199109,199110,199111,199112, \
   199201,199202,199203,199204,199205,199206,199207,199208,199209,199210,199211,199212, \
   199301,199302,199303,199304,199305,199306,199307,199308,199309,199310 /)
H2SO4_Tg = (/ 0.01,5.033,15.092,16.343,21.41,20.418,19.915,17.181,17.363,16.869, \
   18.442,16.526,16.515,18.678,16.025,15.471,15.001,16.729,12.849,13.043,13.008, \
   10.685,9,6.175,5.625,6.674,5.5,4.194,3.722,1.87 /)

caseid=getenv("caseid")

Rearth  = 6.37122e6               ; Earth's radius (m)
PI      = 4.0 * atan(1.0)         ; pi

SAearth=4.e0*PI*Rearth*Rearth ; Earth surface area in m2 (CAM)
;print("SAearth="+SAearth+" m2")
cfac = SAearth/1.e6 ; (kgS/m2 -> GgS)
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
  print("opening "+gmdir+filename)
  a=addfile(gmdir+filename,"r")
  date=a->date
  date=date/100
  date=where(mod(date,100).eq.1,date-89,date-1)
  TMso4_strat   = a->TMso4_strat
  printVarSummary(TMso4_strat)
  time=TMso4_strat&time
  year1=(time(0)-31)/365
  nMonths=dimsizes(TMso4_strat)
  nYears=nMonths/12
  print("nYears="+nYears)
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
    print("=====================================")
    print("Annual Stratospheric burdens (Gg S):")
;    print("OCS   DMS      SO2  H2SO4 sulfate")
;    print(sprintf("%5.1f",TMOCS_strat)+" "+sprintf("%5.1e",TMDMS_strat)+" "+sprintf("%5.2f",TMSO2_strat)+" "+\
;  	  sprintf("%5.3f",TMH2SO4_strat)+" "+sprintf("%5.1f",TMso4_strat))
    print(" Year   OCS   DMS      SO2   H2SO4 sulfate")
    do i = 0,nYears-1
      t1=12*i
      t2=t1+11
      year = date(t1)/100
      print(" "+year +"   " + sprintf("%5.1f",avg(TMOCS_strat(t1:t2)))+" "+\
            sprintf("%5.1e",avg(TMDMS_strat(t1:t2)))+" "+sprintf("%6.2f",avg(TMSO2_strat(t1:t2)))+" "+\
  	    sprintf("%5.3f",avg(TMH2SO4_strat(t1:t2)))+" "+sprintf("%6.1f",avg(TMso4_strat(t1:t2))))
    end do
    
    print("=====================================")
    print("Annual Tropospheric burdens (Gg S):")
;    print("OCS    DMS    SO2  H2SO4 sulfate")
;    print(sprintf("%5.1f",TMOCS_trop)+" "+sprintf("%5.1f",TMDMS_trop)+" "+sprintf("%5.2f",TMSO2_trop)+" "+\
;  	  sprintf("%5.3f",TMH2SO4_trop)+" "+sprintf("%5.1f",TMso4_trop))
    print(" Year   OCS     DMS   SO2    H2SO4  sulfate")
    do i = 0,nYears-1
      t1=12*i
      t2=t1+11
      year = date(t1)/100
      print(" "+year +"   " + sprintf("%6.1f",avg(TMOCS_trop(t1:t2)))+" "+\
            sprintf("%5.1f",avg(TMDMS_trop(t1:t2)))+" "+sprintf("%7.2f",avg(TMSO2_trop(t1:t2)))+" "+\
  	    sprintf("%6.3f",avg(TMH2SO4_trop(t1:t2)))+" "+sprintf("%6.1f",avg(TMso4_trop(t1:t2))))
    end do
    
    b=addfile("/glade/p/work/mmills/inputdata/sulf/CESM_1949_2100_sad_V2_gw_TROPZ_so4strat.nc","r")
    dateb=b->date
    dateb=dateb/100
    indb=ind(dateb.eq.date(0))
;    print("indb="+indb+" date="+dateb(indb))
    CCMI_so4_strat_tot=b->so4_strat_tot
    
    print("=====================================")
    print("Monthly Stratospheric burdens (Gg S):")
    print(" Month    OCS   DMS      SO2     H2SO4  SO4    CCMI Baran/3 Baran/4")
    do i = 0,nMonths-1
      BaranInd=ind(BaranMonth.eq.date(i))
      if (ismissing(BaranInd)) then
        Baran=-1.0
        Baran4=-1.0
      else
        Baran=H2SO4_Tg(BaranInd)*1000.0*32.064999/98.078479
        Baran4=H2SO4_Tg(BaranInd)*1000.0*.75*32.064999/98.078479
      end if
      print(" "+date(i) +"   " + sprintf("%5.1f",TMOCS_strat(i))+" "+\
            sprintf("%5.1e",TMDMS_strat(i))+" "+sprintf("%7.2f",TMSO2_strat(i))+" "+\
  	    sprintf("%6.3f",TMH2SO4_strat(i))+" "+sprintf("%6.1f",TMso4_strat(i))+" "+\
            sprintf("%6.1f",CCMI_so4_strat_tot(indb+i))+" "+sprintf("%6.1f",Baran)+" "+\
            sprintf("%6.1f",Baran4))
    end do
  else
    TMso4_strat=TMso4_strat*MWS/MWH2SO4
    print((/caseid/))
    print("time(0)="+time(0)+" "+TMso4_strat&time@units)
    print( "Stratospheric sulfate burdens (Gg S):" )
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
