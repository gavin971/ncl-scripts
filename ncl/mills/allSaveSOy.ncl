; Constants
;
; NOTE: These are from CAM, shr_const_mod.F90.
BOLTZ   = 1.38065e-23             ; Boltzmann's constant ~ J/K/molecule
AVOGAD  = 6.02214e26              ; Avogadro's number ~ molecules/kmole
RGAS    = AVOGAD*BOLTZ            ; Universal gas constant ~ J/K/kmole
MWDAIR  = 28.966                  ; molecular weight dry air ~ kg/kmole
RDAIR   = RGAS/MWDAIR             ; Dry air gas constant ~ J/K/kg

load "ncl-scripts/getHistoryPath.ncl"

; Calculate parameters of the air based upon the data in the history file.
;
undef("setupair")
procedure setupair(historyFile:file, P, arho)

local hyam, hybm, P0, PS, Q, T, arho

begin
  hyam 	= historyFile->hyam
  hybm 	= historyFile->hybm
  P0 	= historyFile->P0
  PS 	= historyFile->PS
;  Q     = historyFile->Q
  T     = historyFile->T
 
  ; Calculate the air pressures at the interfaces.
  P = pres_hybrid_ccm(PS(:,:,:), P0, hyam(:) ,hybm(:))
  
  ; Calculate the density assuming an ideal gas.
  arho(:,:,:,:)  = P(:,:,:,:) / (RDAIR * T(:,:,:,:))
end

undef("saveSOy")
procedure saveSOy(historyRun:string,pattern:string,subdir:string)

;  inputs: historyRun, i.e. "coag04"
;          pattern, i.e. "nopur2.cam2.h0.1995"

local historyNames,rmin,cpi,rmassmin,writefile,hrc,nhrc,fc

begin
  historyPath = getHistoryPath(historyRun) + subdir
  historyNames = systemfunc("ls " + historyPath + " | grep " + pattern)
  print(historyNames)
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)

  nFiles=dimsizes(historyNames)
  do ifile = 0, nFiles-1
    print("")
    print("    " + ifile + " :  " + (/ historyNames(ifile) /))
    
;    writefile=addfile(historyNames(ifile),"w")
    ; copy history file:
;    allfilename=historyPath+"all."+historyNames(ifile)
    allfilename=historyPath+historyNames(ifile)
    print(allfilename)
;    out = systemfunc("cp " + historyPath + historyNames(ifile) + " " + allfilename )
;    print("file copied")
    writefile=addfile(allfilename,"w")
    
    soy              = writefile->OCS
    mixsulvmr 	     = writefile->MIXSULMD
    dussulmd 	     = writefile->DUSSULMD
    mixsulvmr        = mixsulvmr - dussulmd ; subtract core mass from sulfate (5/19/2009)
    
    if (isfilevar(writefile,"RHOAIR")) then
      arho           = writefile->RHOAIR
      writerho       = False
    else
    ; Compute the air mass and density.
    ;
    ; NOTE: Convert arho and amass are in cgs units.
      P                = soy ; get metadata and array structure from OCS
      P@units          = "Pa"
      P@long_name      = "Pressure"

      arho           = soy ; get metadata and array structure from OCS
      arho@units     = "g/cm3"
      arho@long_name = "air density"

      print("calling setupair")
      setupair(writefile, P, arho) ; calculate pressure and density
      arho = arho(:,:,:,:) * 1e3 / 1e6
      writerho       = True
    end if

    print("      Calculating the mixed sulfate vmr ...")
; (mol/mol)=(g/cm3)/(g/cm3)*(g/mol Air)/(g/mol H2SO4)
;   (28.964 g/mol Air)/(98.0734 g/mol H2SO4) =0.2953
    mixsulvmr=mixsulvmr(:,:,:,:)/arho(:,:,:,:)*0.2953
    mixsulvmr@units      = "mol/mol"
    mixsulvmr@long_name  = "pure sulfates equivalent vmr"

    print("      Calculating the SOy vmr ...")
;    soy      	= writefile->OCS    
    s 	             = writefile->S    
    so 	             = writefile->SO
    so2 	     = writefile->SO2
    so3 	     = writefile->SO3
    hso3             = writefile->HSO3   
    h2so4            = writefile->H2SO4    
    soy=soy+s+so+so2+so3+hso3+h2so4+mixsulvmr
    soy@units        = "mol/mol"
    soy@long_name    = "total sulfur"
    
    if (writerho) then
      writefile->RHOAIR=arho
      writefile->PRESSURE=P
    end if
    
    writefile->MIXSULVMR=mixsulvmr
    
    if (isfilevar(writefile,"PURSULMD")) then
      if (isfilevar(writefile,"PURSULVMR")) then
        pursulvmr=writefile->PURSULVMR
      else
        print("      Calculating the pure sulfate vmr ...")
        pursulvmr 	= writefile->PURSULMD    
        pursulvmr=pursulvmr(:,:,:,:)/arho(:,:,:,:)*0.2953
        pursulvmr@units      = "mol/mol"
        pursulvmr@long_name  = "pure sulfates equivalent vmr"
        writefile->PURSULVMR=pursulvmr
      end if  
      soy=soy+pursulvmr
    end if

    writefile->SOY=soy
  end do
end

;saveSOy("geo2","geo2.cam2.h3.1996-1[12]")
;saveSOy("geo2","geo2.cam2.h3.1997")
;saveSOy("nongeo2","nongeo2.cam2.h3.1997")
;saveSOy("expuT125E01","expuT125E01.cam2.h3.1997-10-17") 
;saveSOy("nogeoPurSpin1","nogeoPurSpin1.cam2.h3.1995-09") 
;saveSOy("nogeoPurSpin1","nogeoPurSpin1.cam2.h3.1995-10") 
;saveSOy("125smin","125smin.cam2.h3") 
;saveSOy("nogeoPurSpin1","all.nogeoPurSpin1.cam2.h3") 
saveSOy("125expuT15b","all.125expuT15b.cam2.h3","daily/")
saveSOy("125expuT15b","all.125expuT15b.cam2.h0","monthly/")
saveSOy("125expuT15","all.125expuT15.cam2.h3","daily/")
saveSOy("125expuT15","all.125expuT15.cam2.h0","monthly/")
saveSOy("expuT125E01","all.expuT125E01.cam2.h3","")
saveSOy("nongeo2","nongeo2.cam2","")
saveSOy("nongeo4","nongeo4.cam2","")
saveSOy("nongeo5","nongeo5.cam2","daily")
saveSOy("nongeo5","nongeo5.cam2","monthly")
saveSOy("geo6","geo6.cam2","daily")
saveSOy("geo6","geo6.cam2","monthly")
