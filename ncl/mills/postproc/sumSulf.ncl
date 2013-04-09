load "../common2.ncl"

;===========================================================
undef("sumSulf")
procedure sumSulf(path:string,filen:string)

begin
  PI = 3.14159
  a=addfile(path+filen,"w")

  mw_so4=115.107340 ; (NH4HSO4) value from chemistry/pp_waccm_mozart_mam3/mo_sim_dat.F90
  mw_air=28.966     ; value from chemistry/modal_aero/modal_aero_newnuc.F90

  if (isfilevar(a,"so4_sum")) then
    print("  so4_sum present in file, skipping sum calculation")
  else
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
  end if
  
  if (isfilevar(a,"PRESSURE")) then
    print("  reading RHOAIR")
    arho=a->RHOAIR
  else
    ; Compute the air mass and density.
    ;
    ; NOTE: Convert arho and amass are in cgs units.
    P               = a->Q ; get metadata and array structure from humidity
    P@units         = "Pa"
    P@long_name     = "Pressure"

    arho            = P ; get metadata and array structure
    arho@units      = "g/cm3"
    arho@long_name  = "air density"

    print("calling setupair")
    setupair(a, P, arho) ; calculate pressure and density
    arho = arho * 1e3 / 1e6

    print("  writing RHOAIR")
    a->RHOAIR=arho
    print("  writing PRESSURE")
    a->PRESSURE=P
    delete(P)
  end if
  
  if (isfilevar(a,"so4_sad")) then
    print("  so4_sad present in file, skipping SAD calculation")
  else
    print("  calculating so4_sad")
    d1=a->dgnw_a01 * 1e6 ; wet diameter of mode 1 in microns
    d2=a->dgnw_a02 * 1e6 ; wet diameter of mode 2 in microns
    d3=a->dgnw_a03 * 1e6 ; wet diameter of mode 3 in microns
    num1=a->num_a1 ; get metadata for so4_sad
    num1=num1 * arho / 1000. ; (#/kg air) * (g air/cm3) / (1000 g/kg)
    num2=a->num_a2 * arho / 1000. ; (#/kg air) * (g air/cm3) / (1000 g/kg)
    num3=a->num_a3 * arho / 1000. ; (#/kg air) * (g air/cm3) / (1000 g/kg)
    so4_sad=num1 ; get metadata
    so4_sad@long_name = "modal aerosol surface area density"  
    so4_sad@units = "ug2/cm2"  
   ; surface area = 4 * PI * (diameter/2)^2 = PI * diameter^2
    so4_sad=PI*(so4_sad*d1*d1+num2*d2*d2+num3*d3*d3)
    print("  writing so4_sad")
    a->so4_sad=so4_sad    
  end if

end ; procedure sumSulf

;===========================================================
; Main driver

caseid = getenv("caseid")

if (ismissing(caseid)) then

  print("set environment variable caseid before running")
  
else

  pathhist="/glade/scratch/mmills/archive/"+caseid+"/atm/hist/"
  path="/glade/scratch/mmills/archive/"+caseid+"/atm/proc/"
  print("mkdir "+path)
  mkdir=systemfunc("mkdir "+path)

  historyNames = systemfunc("ls " + pathhist + " | grep h0")

  nfiles=dimsizes(historyNames)
  do i=0,nfiles-1
    j=i+1
    print("file "+j+" of "+nfiles+":"+historyNames(i))
    if(isfilepresent(path+historyNames(i))) then
      print("  file present in proc")
    else
      print("cp "+pathhist+historyNames(i)+" "+path)
      cp=systemfunc("cp "+pathhist+historyNames(i)+" "+path)
      delete(cp)
    end if
    sumSulf(path,historyNames(i))
  end do
  
end if
