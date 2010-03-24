; Meteoritic Dust - Minimum Temperatures

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

load "ncl-scripts/core.ncl"


  files = systemfunc("ls " + " run/mdust/*cam2.h0.*.nc")
  nFiles = dimsizes(files)
  
  if (nFiles .le. 0) then
    print("ERROR: No files found for ... " + histPath)
  else 
    
    ; Open the last history files
    print("  Opening file ...")
    print("    " + (/ files(:nFiles-1) /))
    
    history = addfiles(files(:nFiles-1), "r")
    ListSetType(history, "cat")
      
    time = history[:]->time  
    T = history[:]->T
    
;;    time = 2025.0 + time/365.0
    time@units     = "days since 2025-1-1"
    time@long_name = "Time"
    
 ;   utc_date = ut_calendar(time, 4)
    
    minT = dim_min(dim_min(dim_min(T)))
    minT@units     = "K"
    minT@long_name = "Temperature (K)"
      
    wks = gsn_open_wks("pdf", "images/min_t")

    nYears = 10
    xbValues = new((/ nYears*4 /), "integer")
    xbLabels = new((/ nYears*4 /), "string")
        
    do i = 0, nYears-1
      xbValues(i*4+0) = 365*i+0 
      xbValues(i*4+1) = 365*i+90 
      xbValues(i*4+2) = 365*i+181 
      xbValues(i*4+3) = 365*i+273
      
      xbLabels(i*4+0) = "J" 
      xbLabels(i*4+1) = "A" 
      xbLabels(i*4+2) = "J:C:" + i 
      xbLabels(i*4+3) = "O" 
    end do
    
    res                 = True
    res@tiMainString    = "Mininum Daily Mean Temperature"
    res@tmXBMode        = "Explicit"
    res@tmXBValues      = xbValues
    res@tmXBLabels      = xbLabels
      
    plot = gsn_csm_xy(wks, time, minT, res)
  end if
  
