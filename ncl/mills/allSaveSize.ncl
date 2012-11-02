load "common2.ncl"
;load "ncl-scripts/getHistoryPath.ncl"

;=========================================

undef("saveSize")
procedure saveSize(historyPath:string,prefix:string,skipbin:integer,sathist:logical)

;  inputs: historyRun, i.e. "coag04"
;          prefix, i.e. "nopur2.cam2.h0.1995"
;          skipbin, missing bin to skip & interpolate (-1 to not skip)

local historyNames,rmin,cpi,rmassmin,writefile,hrc,nhrc,fc

begin

;  historyPath = getHistoryPath(historyRun)
  
;  historyPath = historyPath + "new/"
  
; Create array of files matching pattern
  historyNames = systemfunc("ls " + historyPath + " | grep " + prefix)
  print(historyNames)

  ; Setup the same bin structure used during the model run.
;  nBins = 38
;  rmin  = 2e-8 ; cm
;  rmrat = 2.0
;  rho   = RHO_DUST  ; g/cm-3
  
  nBins = 30
  rmin  = 3.43230298e-8 ; cm
  rmrat = 2.4
  rho   = 1.923  ; g/cm-3
  
  rmass = new(nBins, float)
  rmass!0 = "bins"
  dustrad = rmass
  dustdr = rmass
  rmassup = rmass
  dm = rmass
  dustvol = rmass
  dustradup = rmass
  dustradlow = rmass
  
  print("")
  print("    Assuming ...")
  print("      nBins    =  " + nBins)
  print("      rmin     =  " + rmin)
  print("      rmrat    =  " + rmrat)
  print("      rho      =  " + rho)
  
  setupbins(nBins, rmin, rmrat, rho, rmass, dustrad, dustdr, rmassup, dm, dustvol, dustradup, dustradlow)
  
  rmass@long_name  = "bin center mass"
  rmassup@long_name = "upper bin boundary mass"
  dustrad@long_name = "dust bin mean (volume-weighted) radius"
  dustradup@long_name = "dust upper bin boundary radius"
  dustdr@long_name = "dust bin radius spacing"
  dustradlow@long_name = "dust lower bin boundary radius"
  dustvol@long_name = "dust bin mean volume"
  dm@long_name = "bin mass spacing"
 
  rmass@units = "g"
  rmassup@units = "g"
  dustrad@units = "cm"
  dustradup@units = "cm"
  dustdr@units = "cm"
  dustradlow@units = "cm"
  dustvol@units = "cm3"
  dm@long_name = "g"
    
  
  nFiles=dimsizes(historyNames)
  do ifile = 0, nFiles-1
    print("")
    print("file " + (ifile+1) + " of " + nFiles + " :")
    print("   " + (/ historyNames(ifile) /))
    
    
    ; copy history file:
;    allfilename=historyPath+historyNames(ifile)
    allfilename=historyPath+historyNames(ifile)
;    print("copying file to "+allfilename)
;    out = systemfunc("cp " + historyPath + historyNames(ifile) + " " + allfilename )
;    print("file moved")
    writefile=addfile(allfilename,"w")
    print("file added")

    ; Find the last time.
    lev = writefile->lev
    ilev = writefile->ilev
    nLevs = dimsizes(lev)
    q     = writefile->Q
    
    if (sathist) then
      ncol=dimsizes(q(:,0))
      nFTimes = ncol
    else
      time = writefile->time
      nFTimes = dimsizes(time)
      lat = writefile->lat
      nLats = dimsizes(lat)
      lon = writefile->lon
      nLons = dimsizes(lon)
    end if

    wtpct = writefile->H2SO4WT
    temp  = writefile->T
        
    if (isfilevar(writefile,"SULDEN")) then
      print("copying SULDEN")
      sulfDens=writefile->SULDEN    
    else
      print("calculating SULDEN")
      sulfDens=wtpct
      sulfDens@long_name="sulfate aerosol mass density"
      sulfDens@units="g/cm3"
      print("Calculating sulfate density ...")
      calcSulfDens(wtpct,temp,sulfDens)
      print("Writing SULDEN ...")
      writefile->SULDEN=sulfDens
    end if
    
    if (.not.(isfilevar(writefile,"RMASS"))) then
      print("writing aerosol properties")
      writefile->RMASS=rmass
      writefile->RMASSUP=rmassup
      writefile->DRYRAD=dustrad
      writefile->DRYRADUP=dustradup
      writefile->DRYDR=dustdr
      writefile->DRYRADLOW=dustradlow
      writefile->DRYVOL=dustvol
      writefile->DM=dm
    end if
    
    if (isfilevar(writefile,"PRESSURE")) then
;      P=writefile->PRESSURE
      arho=writefile->RHOAIR
    else
    
    ; Compute the air mass and density.
    ;
    ; NOTE: Convert arho and amass are in cgs units.
      P               = q ; get metadata and array structure from humidity
      P@units         = "Pa"
      P@long_name     = "Pressure"

      arho            = q ; get metadata and array structure from humidity
      arho@units      = "g/cm3"
      arho@long_name  = "air density"

      print("calling setupair")
      setupair(writefile, P, arho) ; calculate pressure and density
      arho = arho * 1e3 / 1e6
      
      writefile->RHOAIR=arho
      writefile->PRESSURE=P
      delete(P)
      
    end if
    
    if (.not.(isfilevar(writefile,"SOY"))) then
      print("reading sulfur constituents ...")
      soy = writefile->OCS    
      s = writefile->S    
      so = writefile->SO
      so2 = writefile->SO2
      so3 = writefile->SO3
      hso3 = writefile->HSO3   
      h2so4 = writefile->H2SO4    

      if (isfilevar(writefile,"MIXSULVMR")) then
        mixsulvmr = writefile->MIXSULVMR
      else
        if (isfilevar(writefile,"MIXSULMD")) then
          print("Calculating the mixed sulfate vmr ...")
          mixsulvmr = writefile->MIXSULMD
          dussulmd = writefile->DUSSULMD
          mixsulvmr = mixsulvmr - dussulmd ; subtract core mass from sulfate (5/19/2009)
    ; (mol/mol)=(g/cm3)/(g/cm3)*(g/mol Air)/(g/mol H2SO4)
    ;   (28.964 g/mol Air)/(98.0734 g/mol H2SO4) =0.2953
          mixsulvmr=mixsulvmr(:,:,:,:)/arho(:,:,:,:)*0.2953
          mixsulvmr@units = "mol/mol"
          mixsulvmr@long_name = "pure sulfates equivalent vmr"
          writefile->MIXSULVMR = mixsulvmr
          delete(dussulmd)
        else
          mixsulvmr=0.0
        end if
      end if
            
      soy=soy+s+so+so2+so3+hso3+h2so4+mixsulvmr
      soy@units      = "mol/mol"
      soy@long_name  = "total sulfur"
      delete(mixsulvmr)
      delete(s)
      delete(so)
      delete(so2)
      delete(so3)
      delete(hso3)
      delete(h2so4)
      
      if (isfilevar(writefile,"PURSULMD")) then
        if (isfilevar(writefile,"PURSULVMR")) then
          pursulvmr = writefile->PURSULVMR
        else
          pursulvmr = writefile->PURSULMD    
          print("Calculating the pure sulfate vmr ...")
          pursulvmr = pursulvmr/arho*0.2953
          pursulvmr@units = "mol/mol"
          pursulvmr@long_name = "pure sulfates equivalent vmr"
          writefile->PURSULVMR = pursulvmr
        end if
        soy=soy+pursulvmr
        delete(pursulvmr)
      end if

      print("writing SOY ...")
      writefile->SOY=soy
      delete(soy)

    end if        

    ; Get the Number density (#/cm-3) for each bin at the 
    ; specified time and location.
    if (isfilevar(writefile,"DUST01") .and. (.not.(isfilevar(writefile,"DUSTNDBIN")))) then
      print("Calculating the dust number density ...")
      if (sathist) then
        dussulndbin = new((/ ncol, nLevs, nBins /), float)
        dussulndbin!0="ncol"
        dussulndbin!1="lev"
        dussulndbin&lev=lev
        dussulndbin!2="bins"   
        dussulndbin@long_name="dust number density per bin"
        dussulndbin@units="/cm3"
      else
        dustndbin = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
        dustndbin!0="time"
        dustndbin&time=time
        dustndbin!1="lev"
        dustndbin&lev=lev
        dustndbin!2="lat"
        dustndbin&lat=lat
        dustndbin!3="lon"
        dustndbin&lon=lon
        dustndbin!4="bins"   
        dustndbin@long_name="dust number density per bin"
        dustndbin@units="/cm3"
      end if
      getN(writefile, "DUST", rmass, arho, dustndbin, skipbin)
     
      print("Writing DUSTNDBIN...")
      writefile->DUSTNDBIN=dustndbin
      delete(dustndbin)
    end if
    
    if (isfilevar(writefile,"PURSUL01") .and. (.not.(isfilevar(writefile,"PURSULNDBIN")))) then
      print("Calculating the pure sulfate number density ...")
      if (sathist) then
        pursulndbin = new((/ ncol, nLevs, nBins /), float)
        pursulndbin!0="ncol"
        pursulndbin!1="lev"
        pursulndbin&lev=lev
        pursulndbin!2="bins"   
        pursulndbin@long_name="pure sulfate number density per bin"
        pursulndbin@units="/cm3"
      else
        pursulndbin = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
        pursulndbin!0="time"
        pursulndbin&time=time
        pursulndbin!1="lev"
        pursulndbin&lev=lev
        pursulndbin!2="lat"
        pursulndbin&lat=lat
        pursulndbin!3="lon"
        pursulndbin&lon=lon
        pursulndbin!4="bins"   
        pursulndbin@long_name="pure sulfate number density per bin"
        pursulndbin@units="/cm3"
      end if
      getN(writefile, "PURSUL", rmass, arho, pursulndbin, skipbin)
      
      print("Writing PURSULNDBIN...")
      writefile->PURSULNDBIN=pursulndbin
      
      delete(pursulndbin)
    end if
 
    if (isfilevar(writefile,"MIXSUL01") .and. (.not.(isfilevar(writefile,"MIXSULNDBIN")))) then
      print("Calculating the mixed sulfate number density ...")
      if (sathist) then
        mixsulndbin = new((/ nFTimes, nLevs, nBins /), float)
        mixsulndbin!0="ncol"
        mixsulndbin!1="lev"
        mixsulndbin&lev=lev
        mixsulndbin!3="bins"   
        mixsulndbin@long_name="mixed sulfate number density per bin"
        mixsulndbin@units="/cm3"
      else
        mixsulndbin = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
        mixsulndbin!0="time"
        mixsulndbin&time=time
        mixsulndbin!1="lev"
        mixsulndbin&lev=lev
        mixsulndbin!2="lat"
        mixsulndbin&lat=lat
        mixsulndbin!3="lon"
        mixsulndbin&lon=lon
        mixsulndbin!4="bins"   
        mixsulndbin@long_name="mixed sulfate number density per bin"
        mixsulndbin@units="/cm3"
      end if
      getN(writefile, "MIXSUL", rmass, arho, mixsulndbin, skipbin)
      
      print("Writing MIXSULNDBIN...")
      writefile->MIXSULNDBIN=mixsulndbin
    end if
 
    if (isfilevar(writefile,"MIXSUL01") .and. ((.not.(isfilevar(writefile,"SULFRAD"))) .or. \
        (.not.(isfilevar(writefile,"CORERAD"))))) then
    
      print("Calculating SULFRAD and CORERAD...")
      
      if (sathist) then
        sulfRad = new((/ ncol, nLevs, nBins /), float)
        sulfRad!0="ncol"
        sulfRad!1="lev"
        sulfRad&lev=lev
        sulfRad!2="bins"
        
        coreRad = new((/ ncol, nLevs, nBins /), float)
        coreRad!0="ncol"
        coreRad!1="lev"
        coreRad&lev=lev
        coreRad!2="bins"   
      else
        sulfRad = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
        sulfRad!0="time"
        sulfRad&time=time
        sulfRad!1="lev"
        sulfRad&lev=lev
        sulfRad!2="lat"
        sulfRad&lat=lat
        sulfRad!3="lon"
        sulfRad&lon=lon
        sulfRad!4="bins"
        
        coreRad = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
        coreRad!0="time"
        coreRad&time=time
        coreRad!1="lev"
        coreRad&lev=lev
        coreRad!2="lat"
        coreRad&lat=lat
        coreRad!3="lon"
        coreRad&lon=lon
        coreRad!4="bins"   
      end if
      sulfRad@long_name="mixed sulfate radius"
      sulfRad@units="cm"
      
      coreRad@long_name="sulfate dust core radius"
      coreRad@units="cm"
      
      calcMixSulRad(writefile, mixsulndbin, wtpct, sulfDens, "DUSSUL", rmass, arho, \
        sulfRad, coreRad, skipbin)
      delete(mixsulndbin)
        
      if (.not.(isfilevar(writefile,"SULFRAD"))) then
        print("saving SULFRAD")
        writefile->SULFRAD=sulfRad
      end if
      delete(sulfRad)
      
      if (.not.(isfilevar(writefile,"CORERAD"))) then
         print("saving CORERAD")
       writefile->CORERAD=coreRad
      end if
      delete(coreRad)
      
    end if
      
    if (isfilevar(writefile,"PURSUL01") .and. (.not.(isfilevar(writefile,"PURERAD")))) then
    
      print("Calculating PURERAD...")
      if (sathist) then
        pureRad = new((/ ncol, nLevs, nBins /), float)
        pureRad!0="ncol"
        pureRad!1="lev"
        pureRad&lev=lev
        pureRad!2="bins"
      else
        pureRad = new((/ nFTimes, nLevs, nLats, nLons, nBins /), float)
        pureRad!0="time"
        pureRad&time=time
        pureRad!1="lev"
        pureRad&lev=lev
        pureRad!2="lat"
        pureRad&lat=lat
        pureRad!3="lon"
        pureRad&lon=lon
        pureRad!4="bins"
      end if
      
      pureRad@long_name="pure sulfate radius"
      pureRad@units="cm"
      
;      if (.not.(isvar(mixsulndbin))) then
;        mixsulndbin = writefile->MIXSULNDBIN
;      end if
      
      calcPurSulRad(writefile, wtpct, sulfDens, rmass, arho, pureRad, skipbin)
      delete(wtpct)
      delete(sulfDens)
      delete(arho)
      
      if (.not.(isfilevar(writefile,"PURERAD"))) then
         print("saving PURERAD")
        writefile->PURERAD=pureRad
      end if
      delete(pureRad)
      
    end if
            
    if (.not.(sathist)) then
      delete(time)
    end if
  end do
end

; Main
; Files located in subdirectory "new" under filepath

;saveSize("125expuT15","125expuT15.cam2.h3.",-1) 
;saveSize("125expuT15","125expuT15.cam2.h3.1997-07-[123]",-1) 
;saveSize("125expuT15","125expuT15.cam2.h3.1997-08",-1) 
;saveSize("125smin","125smin.cam2",-1) 

;saveSize("nogeoPurSpin1","nogeoPurSpin1.cam2",-1) 
;saveSize("nogeo400SO2","nogeo400SO2.cam2",-1) 
;saveSize("125smin400","125smin400.cam2",-1) 
;saveSize("125sminPtbo2","125sminPtbo2.cam2",-1) 

;saveSize("125tabaz22","125tabaz22.cam2",-1) 

;historyPath="/Volumes/Data/Models/cesm/case/f40.Nabro.2deg.wcm.sd.carma.sulf.009/hs/"
historyPath="/glade/proj2/acd/waccm/mmills/runs/f40.Nabro.2deg.wcm.sd.carma.sulf.009c/hs/"
saveSize(historyPath,"f40.Nabro.2deg.wcm.sd.carma.sulf.009.cam.hs.2011-06-12-00000.nc",-1,True)
