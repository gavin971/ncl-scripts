load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

  caseid="f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.QBOon.001"
  sdpath="/glade/scratch/mmills/archive/"+caseid+"/atm/proc/h0zm/"

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
  
  ccmifile="/glade/p/work/mmills/inputdata/CESM_1991_1995_sad_V2_sum_gw_SO2_TROPZ.nc"
  
  a=addfile(ccmifile,"w")
  lat=a->lat
  nLats=dimsizes(lat)
  mon=a->mon
  nTimes=dimsizes(mon)
  aod532=a->aod532
  e532=a->e532
  
  tropz2=linint1_Wrap(tropz&lat,tropz,False,lat,0)
  printVarSummary(tropz2)
  
  tropz3=aod532 ; get metadata
  tropz3@long_name=tropz@long_name
  tropz3@units=tropz@units
  tropz3@_FillValue=tropz@_FillValue
  delete(tropz3@sum_op_ncl)
  
  tropz3(:,:)=tropz3@_FillValue
  tropz3(0:11,:)=tropz2
  tropz3(12:23,:)=tropz2
  tropz3(24:35,:)=tropz2
  tropz3(36:47,:)=tropz2
  tropz3(48:59,:)=tropz2
  tropz3(60,:)=tropz2(0,:)
  printVarSummary(tropz3)
  
  a->TROP_Z=tropz3
  
  aod532strat=aod532 ; get metadata
  aod532strat@long_name = "Stratospheric aerosol optical depth @ 532 nm (vertical integration of e532 for z>=TROP_Z(lat))"
  aod532strat(:,:)=aod532strat@_FillValue
  
  do i = 0,nTimes-1
    print("month "+i)
    do l = 0,nLats-1
      tz=tropz3(i,l)
      aod532strat(i,l)=dim_sum(e532(i, {tz:}, l))
    end do
  end do
  aod532strat(:,:)=aod532strat(:,:)*0.5 ; 0.5 km altitude intervals
  printVarSummary(aod532strat)
  a->aod532strat=aod532strat
  
