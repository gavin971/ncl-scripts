
procedure sumSulf(path:string,filen:string)

begin
  print(filen)
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

  print("  calculating sulfgas")
  sulfgas=a->SO2
  sulfgas=sulfgas+a->H2SO4
  sulfgas=sulfgas+a->DMS
  sulfgas@long_name = "SO2+H2SO4+DMS concentration"
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
caseid="b.e11.B55TRW5CN.f19_g16.MJMtest.005_16-18km_noSO2"

pathhist="/glade/scratch/mmills/archive/"+caseid+"/atm/hist/"
path="/glade/scratch/mmills/archive/"+caseid+"/atm/proc/"
print("mkdir "+path)
mkdir=systemfunc("mkdir "+path)
print("cp "+pathhist+"*h0* "+path)
cp=systemfunc("cp "+pathhist+"*h0* "+path)

historyNames = systemfunc("ls " + path)
print(historyNames)

nfiles=dimsizes(historyNames)
do i=0,nfiles-1
  sumSulf(path,historyNames(i))
end do
