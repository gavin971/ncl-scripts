
procedure sumSulf(path:string,filen:string)

begin
;  print(filen)
  a=addfile(path+filen,"w")

  mw_so4=115.107340 ; (NH4HSO4) value from chemistry/pp_waccm_mozart_mam3/mo_sim_dat.F90
  mw_air=28.966     ; value from chemistry/modal_aero/modal_aero_newnuc.F90

  print("  calculating so4_sum")
  so4_sum=a->so4_a1
  so4_sum=so4_sum+a->so4_a2
  so4_sum=so4_sum+a->so4_a3
  so4_sum=so4_sum*mw_air/mw_so4 ; convert kg/kg to mol/mol
  so4_sum@long_name = "so4 concentration"
  so4_sum@units = "mol/mol"

  print("  writing so4_sum")
  a->so4_sum=so4_sum

  print("  calculating sulfgas & sulfvolc")
  sulfgas=a->SO2
  sulfgas=sulfgas+a->H2SO4
  sulfgas=sulfgas+a->SO3
  sulfgas=sulfgas+a->HSO3
  
  sulfvolc=sulfgas ; get metadata
  sulfvolc=sulfvolc+so4_sum
  sulfvolc@long_name = "SO2+SO3+HSO3+H2SO4+so4 concentration"
  print("  writing sulfvolc")
  a->sulfvolc=sulfvolc
  
  sulfgas=sulfgas+a->DMS
  sulfgas=sulfgas+a->S
  sulfgas=sulfgas+a->SO
  sulfgas=sulfgas+a->OCS
  sulfgas@long_name = "sulfur-bearing gas concentration"
  print("  writing sulfgas")
  a->sulfgas=sulfgas

  print("  calculating sulftot")
  sulftot=sulfgas ; get metadata
  sulftot=sulftot+so4_sum
  sulftot@long_name = "sulfgas+so4_sum concentration"  
  print("  writing sulftot")
  a->sulftot=sulftot

end

;caseid="b.e11.B55TRW5CN.f19_g16.MJMtest.004"
;caseid="b.e11.B55TRW5CN.f19_g16.MJMtest.005_16-18km_noSO2"
caseid = getenv("caseid")

if (ismissing(caseid)) then

  print("set environment variable caseid before running")
  
else

  pathhist="/glade/scratch/mmills/archive/"+caseid+"/atm/hist/"
  path="/glade/scratch/mmills/archive/"+caseid+"/atm/proc/"
  print("mkdir "+path)
  mkdir=systemfunc("mkdir "+path)
;  print("cp "+pathhist+"*h0* "+path)
;  cp=systemfunc("cp "+pathhist+"*h0* "+path)

  historyNames = systemfunc("ls " + pathhist + " | grep h0")
;  print(historyNames)

  nfiles=dimsizes(historyNames)
  do i=0,nfiles-1
    j=i+1
    print("file "+j+" of "+nfiles+":"+historyNames(i))
    if(isfilepresent(path+historyNames(i))) then
      print("  file present in proc, skipping")
    else
      print("cp "+pathhist+historyNames(i)+" "+path)
      cp=systemfunc("cp "+pathhist+historyNames(i)+" "+path)
      sumSulf(path,historyNames(i))
      delete(cp)
    end if
  end do
  
end if
