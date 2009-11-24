load "tests/getHistoryPath.ncl"
load "tests/common2.ncl"

undef("analysisSulfLimit")
procedure analysisSulfLimit(historyRun:string,pattern:string)

;  inputs: historyRun, i.e. "coag04"
;          year, i.e. "1995"

local historyNames, rmin, cpi, rmassmin, readfile, histfile, \
   histfile, time, nTimes, lev, nLevs, lat, nLats, lon, nLons, \
   dustndbin, dustrad, dust05a, dust10a, dust13a, dust30a, \
   pursulndbin, pursulrad, pursul05a, pursul10a, pursul13a, pursul30a, \   
   mixsulndbin, mixsulrad, mixsul05a, mixsul10a, mixsul13a, mixsul30a

begin

  AVG     = 6.02252e+23             ; Avogadro's number (CARMA) ~ molecules/mole
  cpi = 4./3.*PI

  ;  Molecular weights [g/mol]
  gwtH2SO4 = 98.078479

; Determine the path to history files
  historyPath = getHistoryPath(historyRun)
  
; Create array of files matching pattern
  pathNames = systemfunc("ls " + historyPath + "*" + pattern + "*")
  fileNames = systemfunc("ls " + historyPath + "| grep " + pattern)
  
  if (.not.(isfilepresent("analysis/"+historyRun))) then
    a = systemfunc("mkdir analysis/"+historyRun)
  end if
  
  ; Setup the same bin structure used during the model run.
;  nBins = 38
;  nBins = 29
;  rmin  = 2e-8 ; cm
  
  nBins = 36
  rmin  = (gwtH2SO4/AVG/RHO_H2SO4/cpi)^(1./3.) ; cm

  rmrat = 2.0
  rho   = RHO_DUST  ; g/cm-3
  
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
    
  
  nFiles=dimsizes(fileNames)
  do ifile = 0, nFiles-1
    print("")
    print("file " + (ifile+1) + " of " + nFiles + " :")
    print("   " + (/ fileNames(ifile) /))

    if (.not.(isfilepresent("analysis/"+historyRun+"/"+fileNames(ifile)))) then
    
      histfile = addfile(pathNames(ifile), "r")
      analysisFile = addfile("analysis/"+historyRun+"/"+fileNames(ifile), "c")

; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
      filedimdef(analysisFile,"time",-1,True)

      ; Find the last time.
      time = histfile->time
      nTimes = dimsizes(time)
      lev = histfile->lev
      nLevs = dimsizes(lev)
      lat = histfile->lat
      nLats = dimsizes(lat)
      lon = histfile->lon
      nLons = dimsizes(lon)

      temp  = histfile->T

    ; Compute the air mass and density.
    ;
    ; NOTE: Convert arho and amass are in cgs units.
      P               = temp ; get metadata and array structure from temperature
      P@units         = "Pa"
      P@long_name     = "Pressure"

      arho            = temp ; get metadata and array structure from temperature
      arho@units      = "g/cm3"
      arho@long_name  = "air density"
      
      delete(temp)

      print("calling setupair")
      setupair(histfile, P, arho) ; calculate pressure and density
      arho = arho(:,:,:,:) * 1e3 / 1e6

      delete(P)
      
      pursulndbin = new((/ nTimes, nLevs, nLats, nLons, nBins /), float)
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
      pursulrad = pursulndbin
	pursulrad@long_name="pure sulfate radius"
	pursulrad@units="cm"
      mixsulrad = pursulndbin
	mixsulrad@long_name="mixed sulfate radius"
	mixsulrad@units="cm"
      mixsulndbin = pursulndbin
	mixsulndbin@long_name="mixed sulfate number density per bin"
	mixsulndbin@units="/cm3"
      dustndbin = pursulndbin
	dustndbin@long_name="dust number density per bin"
	dustndbin@units="/cm3"

      print("pursul")
      getRad(histfile, "PURSUL", pursulrad, nBins)
      getN(histfile, "PURSUL", rmass, arho, pursulndbin)    

      print("mixsul")
      getRad(histfile, "MIXSUL", mixsulrad, nBins)
      getN(histfile, "MIXSUL", rmass, arho, mixsulndbin)

      print("dust")
      getN(histfile, "DUST", rmass, arho, dustndbin)

      ; get metadata from history
      mixsul05a = histfile->MIXSULND
      mixsul05a(:, :, :, :) = 0.0
      mixsul10a = mixsul05a
      mixsul13a = mixsul05a
      mixsul30a = mixsul05a

      dust05a = mixsul05a
      dust10a = mixsul05a
      dust13a = mixsul05a
      dust30a = mixsul05a

      pursul05a = mixsul05a
      pursul10a = mixsul05a
      pursul13a = mixsul05a
      pursul30a = mixsul05a
      pursulnd = mixsul05a

      mixsul05a@long_name="mixed sulfates > 0.5 nm number density"
      mixsul10a@long_name="mixed sulfates > 1.0 nm number density"
      mixsul13a@long_name="mixed sulfates > 1.3 nm number density"
      mixsul30a@long_name="mixed sulfates > 3.0 nm number density"

      dust05a@long_name="dust > 0.5 nm number density"
      dust10a@long_name="dust > 1.0 nm number density"
      dust13a@long_name="dust > 1.3 nm number density"
      dust30a@long_name="dust > 3.0 nm number density"

      pursul05a@long_name="pure sulfates > 0.5 nm number density"
      pursul10a@long_name="pure sulfates > 1.0 nm number density"
      pursul13a@long_name="pure sulfates > 1.3 nm number density"
      pursul30a@long_name="pure sulfates > 3.0 nm number density"

      ; tally particle concentrations greater than .5, 1.0, 1.3. and 3 nm 
      ; NOTE: this tally could increase due to partial bin counts, not included yet...

      mixThresh=False  ; Flag for when smallest sulfates are greater than largest threshold
      pureThresh=False  ; Flag for when smallest sulfates are greater than largest threshold

      sec2=stringtointeger(systemfunc("date +%s"))
      print("looping bins 0 to "+(nBins-1)+"...")
      do iBin = 0, nBins-1

	if (dustrad(iBin).ge.5e-8) then
          dust05a(:, :, :, :) = dust05a(:, :, :, :) + dustndbin(:, :, :, :, iBin)
	end if
	if (dustrad(iBin).ge.1e-7) then
          dust10a(:, :, :, :) = dust10a(:, :, :, :) + dustndbin(:, :, :, :, iBin)
	end if
	if (dustrad(iBin).ge.1.3e-7) then
          dust13a(:, :, :, :) = dust13a(:, :, :, :) + dustndbin(:, :, :, :, iBin)
	end if
	if (dustrad(iBin).ge.3e-7) then
          dust30a(:, :, :, :) = dust30a(:, :, :, :) + dustndbin(:, :, :, :, iBin)
	end if


	if (.not. mixThresh) then
          if (min(mixsulrad(:, :, :, :, iBin)).ge.3.e-7) then
            mixThresh=True
          end if
	end if

	if (mixThresh) then
          mixsul05a(:, :, :, :) = mixsul05a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
          mixsul10a(:, :, :, :) = mixsul10a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
          mixsul13a(:, :, :, :) = mixsul13a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
          mixsul30a(:, :, :, :) = mixsul30a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)  
	else     
          mixsul05a(:, :, :, :) = mixsul05a(:, :, :, :) + where(mixsulrad(:, :, :, :, iBin).ge.5e-8,   mixsulndbin(:, :, :, :, iBin), 0.0)
          mixsul10a(:, :, :, :) = mixsul10a(:, :, :, :) + where(mixsulrad(:, :, :, :, iBin).ge.1e-7,   mixsulndbin(:, :, :, :, iBin), 0.0)
          mixsul13a(:, :, :, :) = mixsul13a(:, :, :, :) + where(mixsulrad(:, :, :, :, iBin).ge.1.3e-7, mixsulndbin(:, :, :, :, iBin), 0.0)
          mixsul30a(:, :, :, :) = mixsul30a(:, :, :, :) + where(mixsulrad(:, :, :, :, iBin).ge.3e-7,   mixsulndbin(:, :, :, :, iBin), 0.0)
	end if


	if (.not. pureThresh) then
          if (min(pursulrad(:, :, :, :, iBin)).ge.3.e-7) then
            pureThresh=True
          end if
	end if

	if (pureThresh) then
          pursul05a(:, :, :, :) = pursul05a(:, :, :, :) + pursulndbin(:, :, :, :, iBin)
          pursul10a(:, :, :, :) = pursul10a(:, :, :, :) + pursulndbin(:, :, :, :, iBin)
          pursul13a(:, :, :, :) = pursul13a(:, :, :, :) + pursulndbin(:, :, :, :, iBin)
          pursul30a(:, :, :, :) = pursul30a(:, :, :, :) + pursulndbin(:, :, :, :, iBin)  
	else
          pursul05a(:, :, :, :) = pursul05a(:, :, :, :) + where(pursulrad(:, :, :, :, iBin).ge.5e-8,   pursulndbin(:, :, :, :, iBin), 0.0)
          pursul10a(:, :, :, :) = pursul10a(:, :, :, :) + where(pursulrad(:, :, :, :, iBin).ge.1e-7,   pursulndbin(:, :, :, :, iBin), 0.0)
          pursul13a(:, :, :, :) = pursul13a(:, :, :, :) + where(pursulrad(:, :, :, :, iBin).ge.1.3e-7, pursulndbin(:, :, :, :, iBin), 0.0)
          pursul30a(:, :, :, :) = pursul30a(:, :, :, :) + where(pursulrad(:, :, :, :, iBin).ge.3e-7,   pursulndbin(:, :, :, :, iBin), 0.0)
	end if

	sec=sec2
	sec2=stringtointeger(systemfunc("date +%s"))
	sec=sec2-sec
	print("bin "+iBin+": "+sec+" sec")
      end do ; iBin

      delete(mixsulndbin)
      delete(mixsulrad)

      delete(dustndbin)

      delete(pursulndbin)
      delete(pursulrad)

      delete(sec)
      delete(sec2)

      print("Calculating zonal means...")

      print("pursul")
      
      pursul05a_zm=dim_avg_Wrap(pursul05a(time|:,lev|:,lat|:,lon|:))
      pursul05a_zm@long_name="Pure sulfate number density (r>0.5 nm), zonal mean"
      
      pursul10a_zm=dim_avg_Wrap(pursul10a(time|:,lev|:,lat|:,lon|:))
      pursul10a_zm@long_name="Pure sulfate number density (r>1 nm), zonal mean"
      
      pursul13a_zm=dim_avg_Wrap(pursul13a(time|:,lev|:,lat|:,lon|:))
      pursul13a_zm@long_name="Pure sulfate number density (r>1.3 nm), zonal mean"
      
      pursul30a_zm=dim_avg_Wrap(pursul30a(time|:,lev|:,lat|:,lon|:))
      pursul30a_zm@long_name="Pure sulfate number density (r>3 nm), zonal mean"
      
      print("Saving pursul")
      analysisFile->pursul05a_zm=pursul05a_zm
      analysisFile->pursul10a_zm=pursul10a_zm
      analysisFile->pursul13a_zm=pursul13a_zm
      analysisFile->pursul30a_zm=pursul30a_zm

      print("mixsul")
      mixsul05a_zm=dim_avg_Wrap(mixsul05a(time|:,lev|:,lat|:,lon|:))
      mixsul05a_zm@long_name="Mixed sulfate number density (r>0.5 nm), zonal mean"
      
      mixsul10a_zm=dim_avg_Wrap(mixsul10a(time|:,lev|:,lat|:,lon|:))
      mixsul10a_zm@long_name="Mixed sulfate number density (r>1 nm), zonal mean"
      
      mixsul13a_zm=dim_avg_Wrap(mixsul13a(time|:,lev|:,lat|:,lon|:))
      mixsul13a_zm@long_name="Mixed sulfate number density (r>1.3 nm), zonal mean"
      
      mixsul30a_zm=dim_avg_Wrap(mixsul30a(time|:,lev|:,lat|:,lon|:))
      mixsul30a_zm@long_name="Mixed sulfate number density (r>3 nm), zonal mean"
      
      print("Saving mixsul")
      analysisFile->mixsul05a_zm=mixsul05a_zm
      analysisFile->mixsul10a_zm=mixsul10a_zm
      analysisFile->mixsul13a_zm=mixsul13a_zm
      analysisFile->mixsul30a_zm=mixsul30a_zm

      print("dust")
      dust05a_zm=dim_avg_Wrap(dust05a(time|:,lev|:,lat|:,lon|:))
      dust10a_zm=dim_avg_Wrap(dust10a(time|:,lev|:,lat|:,lon|:))
      dust13a_zm=dim_avg_Wrap(dust13a(time|:,lev|:,lat|:,lon|:))
      dust30a_zm=dim_avg_Wrap(dust30a(time|:,lev|:,lat|:,lon|:))
      print("Saving dust")
      analysisFile->dust05a_zm=dust05a_zm
      analysisFile->dust10a_zm=dust10a_zm
      analysisFile->dust13a_zm=dust13a_zm
      analysisFile->dust30a_zm=dust30a_zm
      
      print("Saving pressure surfaces...")
      print("0.007 hPa")
      print("pursul")
      iLev=0.007
      pursul05a_007hPa=pursul05a(:,{iLev},:,:)
      pursul05a_007hPa@long_name="Pure sulfate number density (r>0.5 nm) at 0.007hPa"      
      analysisFile->pursul05a_007hPa=pursul05a_007hPa
      
      pursul10a_007hPa=pursul10a(:,{iLev},:,:)
      pursul10a_007hPa@long_name="Pure sulfate number density (r>1 nm) at 0.007hPa"      
      analysisFile->pursul10a_007hPa=pursul10a_007hPa
      
      pursul13a_007hPa=pursul13a(:,{iLev},:,:)
      pursul13a_007hPa@long_name="Pure sulfate number density (r>1.3 nm) at 0.007hPa"      
      analysisFile->pursul13a_007hPa=pursul13a_007hPa
      
      pursul30a_007hPa=pursul30a(:,{iLev},:,:)
      pursul30a_007hPa@long_name="Pure sulfate number density (r>3 nm) at 0.007hPa"      
      analysisFile->pursul30a_007hPa=pursul30a_007hPa

      print("mixsul")
      mixsul05a_007hPa=mixsul05a(:,{iLev},:,:)
      mixsul05a_007hPa@long_name="Mixed sulfate number density (r>0.5 nm) at 0.007hPa"      
      analysisFile->mixsul05a_007hPa=mixsul05a_007hPa
      
      mixsul10a_007hPa=mixsul10a(:,{iLev},:,:)
      mixsul10a_007hPa@long_name="Mixed sulfate number density (r>1 nm) at 0.007hPa"      
      analysisFile->mixsul10a_007hPa=mixsul10a_007hPa
      
      mixsul13a_007hPa=mixsul13a(:,{iLev},:,:)
      mixsul13a_007hPa@long_name="Mixed sulfate number density (r>1.3 nm) at 0.007hPa"      
      analysisFile->mixsul13a_007hPa=mixsul13a_007hPa
      
      mixsul30a_007hPa=mixsul30a(:,{iLev},:,:)
      mixsul30a_007hPa@long_name="Mixed sulfate number density (r>3 nm) at 0.007hPa"      
      analysisFile->mixsul30a_007hPa=mixsul30a_007hPa

      print("dust")
      dust05a_007hPa=dust05a(:,{iLev},:,:)
      dust05a_007hPa@long_name="Dust number density (r>0.5 nm) at 0.007hPa"      
      analysisFile->dust05a_007hPa=dust05a_007hPa
      
      dust10a_007hPa=dust10a(:,{iLev},:,:)
      dust10a_007hPa@long_name="Dust number density (r>1 nm) at 0.007hPa"      
      analysisFile->dust10a_007hPa=dust10a_007hPa
      
      dust13a_007hPa=dust13a(:,{iLev},:,:)
      dust13a_007hPa@long_name="Dust number density (r>1.3 nm) at 0.007hPa"      
      analysisFile->dust13a_007hPa=dust13a_007hPa
      
      dust30a_007hPa=dust30a(:,{iLev},:,:)
      dust30a_007hPa@long_name="Dust number density (r>3 nm) at 0.007hPa"      
      analysisFile->dust30a_007hPa=dust30a_007hPa
      
      print("cleaning up...")

      delete(histfile)

      delete(time)
      delete(nTimes)
      delete(lev)
      delete(nLevs)
      delete(lat)
      delete(nLats)
      delete(lon)
      delete(nLons)
    
    end if
    print("returning...")
   
  end do ; ifile
  
  delete(dustrad)
  
end


; Main

;analysisSulfLimit("125kelvH2Og","h3")
analysisSulfLimit("125bin36c","h3")
