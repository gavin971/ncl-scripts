; PMC Temperatures
;
; Comparisson to plots from Hunten1980jas

; Constants
;
; NOTE: These are from CAM, shr_const_mod.F90.
PI      = 3.14159265358979323846  ; pi
G       = 9.80616                 ; acceleration of gravity ~ m/s^2
BOLTZ   = 1.38065e-23             ; Boltzmann's constant ~ J/K/molecule
AVOGAD  = 6.02214e26              ; Avogadro's number ~ molecules/kmole
RGAS    = AVOGAD*BOLTZ            ; Universal gas constant ~ J/K/kmole
MWDAIR  = 28.966                  ; molecular weight dry air ~ kg/kmole
RDAIR   = RGAS/MWDAIR             ; Dry air gas constant ~ J/K/kg
MWH2O   = 18.000                  ; molecular weight of water vapor
PTOP    = 100.0                   ; pressure at the top of the model (non-zero) ~ Pa
REARTH  = 6.37122e6               ; Earth's radius (m)

load "tests/core.ncl"


undef("getT_1999")
procedure getT_1999(WACCM_RUN:string)

local files, nFiles, histories

begin

  ; Get date from April thorugh October.
  histPath = "run/" + WACCM_RUN + "/" + WACCM_RUN + "*cam2.h3.1997-0[456789]*.nc"
  files = systemfunc("ls " + histPath)
  nFiles = dimsizes(files)
  
  if (nFiles .le. 0) then
    print("ERROR: No files found for ... " + histPath)
  else 
  
    ; Get the total number of times.
    print("")
    print("  Opening Files ...")
    print("    " + files)
    histories = addfiles(files, "r")
    ListSetType(histories, "cat")
    
    P0    = histories[0]->P0
    hyam  = histories[0]->hyam
    hybm  = histories[0]->hybm
    
    
    time  = addfiles_GetVar(histories, files, "time")
    T     = addfiles_GetVar(histories, files, "T")
    omega = addfiles_GetVar(histories, files, "OMEGA")
    PS    = addfiles_GetVar(histories, files, "PS")
    Q     = addfiles_GetVar(histories, files, "Q")
    U     = addfiles_GetVar(histories, files, "U")
    V     = addfiles_GetVar(histories, files, "V")
    Z3    = addfiles_GetVar(histories, files, "Z3")
    
    ; Convert omega to w
    print("")
    print("  Calculating P and W ...")
    P = pres_hybrid_ccm(PS, P0, hyam, hybm)
    rho = P / (RGAS * T)
    W = omega
    W = (/ - omega / (rho * G) /) 

    tlat = T(:, :, {70.0}, :)
    qlat = Q(:, :, {70.0}, :) * (MWDAIR / MWH2O) * 1e6
    ulat = U(:, :, {70.0}, :)
    vlat = V(:, :, {70.0}, :)
    wlat = W(:, :, {70.0}, :)
    zlat = Z3(:, :, {70.0}, :) / 1e3
 
    nTimes = dimsizes(time)
    nLons = dimsizes(T&lon)
    
    ; Create an output Z structure
    zs = fspan(35.0, 100., 130)
    zs@long_name = "Altitude"
    zs@units = "km"
    nZs = dimsizes(zs)

    ;Define the output fields
    tbyz2 = new((/ nTimes, nZs, nLons /), "float")
    qbyz2 = tbyz2
    ubyz2 = tbyz2
    vbyz2 = tbyz2
    wbyz2 = tbyz2
    
    print("")
    print("  Interpolating to z at 70N ...")
    do i = 0, nTimes-1
      do j = 0, nLons-1
        tbyz2(i, :, j) = linint1(zlat(i, ::-1, j), tlat(i, ::-1, j), False, zs, 0)
        qbyz2(i, :, j) = linint1(zlat(i, ::-1, j), qlat(i, ::-1, j), False, zs, 0)
        ubyz2(i, :, j) = linint1(zlat(i, ::-1, j), ulat(i, ::-1, j), False, zs, 0)
        vbyz2(i, :, j) = linint1(zlat(i, ::-1, j), vlat(i, ::-1, j), False, zs, 0)
        wbyz2(i, :, j) = linint1(zlat(i, ::-1, j), wlat(i, ::-1, j), False, zs, 0)
      end do
    end do
    
    tbyz = dim_avg(tbyz2)
    tbyz!0 = "time"
    tbyz&time = time
    tbyz!1 = "z"
    tbyz&z = zs
    tbyz@units = "K"
    tbyz@long_name = "Temperature"
    
    qbyz = tbyz
    qbyz = (/ dim_avg(qbyz2) /)
    qbyz@units = "ppmv"
    qbyz@long_name = "Water Vapor"
   
    ubyz = tbyz
    ubyz = (/ dim_avg(ubyz2) /)
    ubyz@units = "m/s"
    ubyz@long_name = "Zonal Wind"
    
    vbyz = tbyz
    vbyz = (/ dim_avg(vbyz2) /)
    vbyz@units = "m/s"
    vbyz@long_name = "Meridional Wind"
    
    wbyz = tbyz
    wbyz = (/ dim_avg(wbyz2) /)
    wbyz@units = "m/s"
    wbyz@long_name = "Vertical Wind"
    
;    print("  " + zs + "  " + wbyz(0, :))
    print("  " + time + "  " + tbyz(:, {82}) + "  " + qbyz(:, {82}) + "  " + wbyz(:, {82}))
       

    tfile = "tlub_1999.nc"
    system("rm -f analysis/" + WACCM_RUN + "/" + tfile) 
    outfile       = addfile("analysis/" + WACCM_RUN + "/" + tfile, "c")
    outfile->time = time
    outfile->z    = zs
    outfile->T    = tbyz
    outfile->Q    = qbyz
    outfile->U    = ubyz
    outfile->V    = vbyz
    outfile->W    = wbyz
  end if    
end


undef("getT_2004")
procedure getT_2004(WACCM_RUN:string)

local files, nFiles, histories

begin

  ; Get date from April thorugh October.
  histPath = "run/" + WACCM_RUN + "/*cam2.h2.*.nc"
  files = systemfunc("ls " + histPath)
  nFiles = dimsizes(files)
  
  ; We want the second year's january and February
  firstFile = 30
  
  if (nFiles .le. 0) then
    print("ERROR: No files found for ... " + histPath)
  else 
  
    ; Get the total number of times.
    histories = addfiles(files(firstFile:), "r")
    ListSetType(histories, "cat")
    
    time = addfiles_GetVar(histories, files(firstFile:), "time")
    T    = addfiles_GetVar(histories, files(firstFile:), "T")
    Z3   = addfiles_GetVar(histories, files(firstFile:), "Z3")

    tlat = T(:, :, {-68.0}, :)
    zlat = Z3(:, :, {-68.0}, :) / 1e3
 
    nTimes = dimsizes(time)
    nLons = dimsizes(T&lon)
    
    ; Create an output Z structure
    zs = fspan(35.0, 100., 130)
    zs@long_name = "Altitude"
    zs@units = "km"
    nZs = dimsizes(zs)

    ;Define the output fields
    tbyz2 = new((/ nTimes, nZs, nLons /), "float")
    
    do i = 0, nTimes-1
      do j = 0, nLons-1
        tbyz2(i, :, j) = linint1(zlat(i, ::-1, j), tlat(i, ::-1, j), False, zs, 0)
      end do
    end do
    
    tbyz = dim_avg(tbyz2)
    tbyz!0 = "time"
    tbyz&time = time
    tbyz!1 = "z"
    tbyz&z = zs
    tbyz@units = "K"
    tbyz@long_name = "Temperature"
      

    tfile = "tlub_2004.nc"
    system("rm analysis/" + WACCM_RUN + "/" + tfile) 
    outfile       = addfile("analysis/" + WACCM_RUN + "/" + tfile, "c")
    outfile->time = time
    outfile->z    = zs
    outfile->T    = tbyz
  end if    
end


; Main
;WACCM_RUN = getenv("WACCM_RUN")

;if (ismissing(WACCM_RUN)) then
;  WACCM_RUN = "results"
;end if
;WACCM_RUN = "expuT125E01"
WACCM_RUN = "ep10rt02"
print("Using WACCM_RUN: " + WACCM_RUN)
print("")

getT_1999(WACCM_RUN)
;getT_2004(WACCM_RUN)

;getT_1999("gw125_e1")
;getT_1999("gw15_e0875")
;getT_1999("gw15_e0875_smax")
;getT_1999("pmc_smin")

;getT_2004("gw15_e0875")
;getT_2004("pmc_smin")
