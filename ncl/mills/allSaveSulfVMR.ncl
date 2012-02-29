; Constants
;
; NOTE: These are from CAM, shr_const_mod.F90.
BOLTZ   = 1.38065e-23             ; Boltzmann's constant ~ J/K/molecule
AVOGAD  = 6.02214e26              ; Avogadro's number ~ molecules/kmole
RGAS    = AVOGAD*BOLTZ            ; Universal gas constant ~ J/K/kmole
MWDAIR  = 28.966                  ; molecular weight dry air ~ kg/kmole
RDAIR   = RGAS/MWDAIR             ; Dry air gas constant ~ J/K/kg

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "ncl-scripts/getHistoryPath.ncl"

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

undef("saveSulfVMR")
procedure saveSulfVMR(historyRun:string,pattern:string)

;  inputs: historyRun, i.e. "coag04"
;          pattern, i.e. "nopur2.cam2.h0.1995"

local historyNames,rmin,cpi,rmassmin,writefile,hrc,nhrc,fc

begin
  host   = systemfunc("hostname")
  if (host.eq."modeling1") then
    rootPath="/data4/mmills/"
  else
    rootPath="$MODEL_PATH/ccsm/case/"
  end if
  
  historyPath = rootPath+historyRun+"/h0/"
  outputPath  = rootPath+historyRun+"/proc/sulvmr/"
  a = systemfunc("mkdir -p "+outputPath)
  
  historyNames = systemfunc("ls " + historyPath + " | grep " + pattern )
  print(historyNames)
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)

  nFiles=dimsizes(historyNames)
  do ifile = 0, nFiles-1
    print("")
    print("    " + ifile + " :  " + (/ historyNames(ifile) /))
    
    writefilename=outputPath+historyNames(ifile)
    
    if (.not.(isfilepresent(writefilename))) then
      readfilename=historyPath+historyNames(ifile)

      readfile =addfile(readfilename,"r")

      pursulvmr 	= readfile->PURSULMD 
        
      if (isfilevar(readfile,"RHOAIR")) then
        arho           = readfile->RHOAIR
        writerho       = False
      else
      ; Compute the air mass and density.
      ;
      ; NOTE: Convert arho and amass are in cgs units.
        P                = pursulvmr ; get metadata and array structure
        P@units          = "Pa"
        P@long_name      = "Pressure"

        arho           = pursulvmr ; get metadata and array structure
        arho@units     = "g/cm3"
        arho@long_name = "air density"

        print("calling setupair")
        setupair(readfile, P, arho) ; calculate pressure and density
        arho = arho(:,:,:,:) * 1e3 / 1e6
        writerho       = True
      end if

      print("creating "+writefilename)
      writefile=addfile(writefilename,"c")
      if (writerho) then
        writefile->RHOAIR=arho
        writefile->PRESSURE=P
      end if

      print("      Calculating the pure sulfate vmr ...")
      pursulvmr=pursulvmr(:,:,:,:)/arho(:,:,:,:)*0.2953
      pursulvmr@units      = "mol/mol"
      pursulvmr@long_name  = "pure sulfates equivalent vmr"
      writefile->PURSULVMR=pursulvmr
      
      if (isfilevar(readfile,"MIXSULMD")) then
        mixsulvmr 	 = readfile->MIXSULMD
        mixsulshellvmr   = mixsulvmr
        dussulmd 	 = readfile->DUSSULMD
        mixsulshellvmr   = mixsulshellvmr - dussulmd ; subtract core mass from sulfate (5/19/2009)

        print("      Calculating the mixed sulfate vmr ...")
    ; (mol/mol)=(g/cm3)/(g/cm3)*(g/mol Air)/(g/mol H2SO4)
    ;   (28.964 g/mol Air)/(98.0734 g/mol H2SO4) =0.2953
        mixsulvmr=mixsulvmr(:,:,:,:)/arho(:,:,:,:)*0.2953
        mixsulshellvmr=mixsulshellvmr(:,:,:,:)/arho(:,:,:,:)*0.2953
        mixsulvmr@units      = "mol/mol"
        mixsulvmr@long_name  = "mixed sulfates equivalent vmr"
        mixsulshellvmr@units      = "mol/mol"
        mixsulshellvmr@long_name  = "mixed sulfates shell equivalent vmr"

        writefile->MIXSULVMR=mixsulvmr
        writefile->MIXSULSHELLVMR=mixsulshellvmr
      end if

    end if
    
  end do
  
end

undef("concatN2OSulfVMR")
procedure concatN2OSulfVMR(historyRun:string,pattern:string)

;  inputs: historyRun, i.e. "coag04"
;          pattern, i.e. "nopur2.cam2.h0.1995"

local historyNames

begin

  host   = systemfunc("hostname")
  if (host.eq."modeling1") then
    rootPath="/data4/mmills/"
  else
    rootPath="$MODEL_PATH/ccsm/case/"
  end if
  
  historyPath = rootPath+historyRun+"/h0/"
  sulfPath    = rootPath+historyRun+"/proc/sulvmr/"
  outputPath  = rootPath+historyRun+"/proc/"
  
  historyNames = historyPath + "/" + systemfunc("ls " + historyPath + " | grep " + pattern )
  print(historyNames)

  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  N2O = addfiles_GetVar(history, historyNames, "N2O")

  historyNames2 = sulfPath + "/" + systemfunc("ls " + sulfPath + " | grep " + pattern)
  print(historyNames2)
  history2 = addfiles(historyNames2, "r")
  ListSetType(history2, "cat")
  PURSULVMR = addfiles_GetVar(history2, historyNames2, "PURSULVMR")
  
  writefile=addfile(outputPath+historyRun+".N2OSulfVMR."+pattern+".nc", "c")
  
  writefile->N2O = N2O
  writefile->PURSULVMR = PURSULVMR
  
  readfile =addfile(historyNames2(0),"r")
  if (isfilevar(readfile,"MIXSULVMR")) then
    MIXSULVMR = addfiles_GetVar(history2, historyNames2, "MIXSULVMR")
    MIXSULSHELLVMR = addfiles_GetVar(history2, historyNames2, "MIXSULSHELLVMR")

    writefile->MIXSULVMR = MIXSULVMR
    writefile->MIXSULSHELLVMR = MIXSULSHELLVMR
  end if
  
end

;saveSulfVMR("f40.2008-2011.2deg.wcm.sd.carma.sulf.002","cam2.h0")
concatN2OSulfVMR("f40.2008-2011.2deg.wcm.sd.carma.sulf.002","cam2.h0")
