load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Constants
;
; NOTE: These are from CAM, shr_const_mod.F90.
MWH2SO4 = 98.0784                 ; molecular weight of H2SO4 species (g/mol)
MWS     = 32.066                  ; molecular weight of atomic sulfur (g/mol)
AVOGAD  = 6.02214e23              ; Avogadro's number ~ molecules/mole
Rearth  = 6.37122e6               ; Earth's radius (m)
PI      = 4.0 * atan(1.0)         ; pi

;  caseid="f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.QBOon.001"
;  sdpath="/glade/scratch/mmills/archive/"+caseid+"/atm/proc/h0zm/"

  caseid="f.e12.FSDW5.f19_f19.1990spinup.001"
  sdpath="/glade/p/cesm/wawg/mmills/run/"+caseid+"/atm/proc/h0zm/"

  print("historyNames = systemfunc(ls "+ sdpath + "*h0zm.1992*))")
  historyNames = systemfunc("ls " + sdpath + "*h0zm.1992*")
  print(historyNames)
  print("history = addfiles(historyNames, r)")
  history = addfiles(historyNames, "r")
  print("ListSetType(history, cat)")
  ListSetType(history, "cat")
  ;  print(historyNames)
  print("historyNames:"+dimsizes(historyNames))

  print ("reading time")
  time    = addfiles_GetVar(history, historyNames, "time")
  
  print ("reading date")
  date    = addfiles_GetVar(history, historyNames, "date")
  print(time+" "+date)
  
  print ("reading TROP_Z")
  tropz    = addfiles_GetVar(history, historyNames, "TROP_Z")
  tropz=tropz/1000
  tropz@units="km"
  printVarSummary(tropz)
  
  printVarSummary(tropz&lat)
    
;  ccmifile="/glade/p/work/mmills/inputdata/CESM_1991_1995_sad_V2_sum_gw_SO2_TROPZ.nc"
  ccmipath="/glade/p/work/mmills/inputdata/sulf/"
  filein ="CESM_1949_2100_sad_V2_c130627.nc"
  fileout="CESM_1949_2100_sad_V2_gw_TROPZ_so4strat.nc"
  command="cp "+ccmipath+filein+" "+ccmipath+fileout
  print((/command/))
  dum=systemfunc(command)
  print(dum)
  
  a=addfile(ccmipath+fileout,"w")
  
  lat=a->lat
  nLats=dimsizes(lat)
  mon=a->mon
  date2=a->date
  nTimes=dimsizes(mon)
  e532=a->e532
  
  gw = lat ; get metadata
  gw = (/ cos(lat*PI/180.) /)
  
  tropz2=linint1_Wrap(tropz&lat,tropz,False,lat,0)
  tropz2!0="month"
  printVarSummary(tropz)
  printVarSummary(tropz2)
  
  aod532=dim_sum_n_Wrap(e532,1)
  aod532@long_name="Aerosol optical depth @ 532 nm (vertical integration of e532)"
  
  a->TROP_Z=tropz2
  a->gw=gw  
  
  so4=a->H2SO4_mass ; molecules/cm3
  factor=1e3*MWS/AVOGAD ; molecules/cm3 -> g/cm3 * 1e6/1e3 -> kg/m3
  so4=so4*factor ; kg/m3
  
  so4_strat=aod532 ; get metadata and array structure
  so4_strat@long_name="Stratospheric sulfate mass"
  so4_strat@units="kg-S/m2"
  so4_strat(:,:)=so4_strat@_FillValue
  
  aod532strat=aod532 ; get metadata
  aod532strat@long_name = "Stratospheric aerosol optical depth @ 532 nm (vertical integration of e532 for z>=TROP_Z(lat))"
  aod532strat(:,:)=aod532strat@_FillValue
  
  do i = 0,nTimes-1
    i1=i+1
    print("month "+i1+" of "+nTimes+": "+date2(i)+" "+mon(i))
    monind=floattoint(mon(i)-1)
    do l = 0,nLats-1
      tz=tropz2(monind,l)
      if (ismissing(tz)) then
        print("lat "+l+": tz="+tz)
      else
        aod532strat(i,l)=dim_sum(e532(i, {tz:}, l))
        so4_strat(i,l)=dim_sum(so4(i, {tz:}, l))
      end if
   end do
  end do
  aod532strat(:,:)=aod532strat(:,:)*0.5 ; 0.5 km altitude intervals
  printVarSummary(aod532strat)
  a->aod532strat=aod532strat
  
  so4_strat(:,:)=so4_strat(:,:)*500.e0 ; 0.5 km altitude intervals
  a->so4_strat=so4_strat

  SAearth=4.e0*PI*Rearth*Rearth ; Earth surface area in m2 (CAM)
  print("SAearth="+SAearth+" m2")
  cfac = SAearth/1.e6 ; (kgS/m2 -> TgS)
  print("cfac="+cfac+" kgS/m2 -> TgS")
  
  so4_strat_tot=dim_avg_wgt_Wrap(so4_strat,gw,0)
  so4_strat_tot=so4_strat_tot*cfac
  so4_strat_tot@long_name="Global stratospheric sulfate mass"
  so4_strat@units="Tg-S"
  a->so4_strat_tot=so4_strat_tot
  
  ; For a nicer looking netCDF, create a "new line" character.
  nl = integertochar(10)  ; newline character
  
  globalAtt             = True
  globalAtt@creation_date     = nl+\
      systemfunc("date")
  globalAtt@creator = nl+\
      "Mike Mills, mmills@ucar.edu"
  globalAtt@history = nl+\
       "Adapted from file "+filein+\
       "Additional variables added by program addCCMITropZ.ncl"
  fileattdef( a, globalAtt )
  
  delete(a)
  
