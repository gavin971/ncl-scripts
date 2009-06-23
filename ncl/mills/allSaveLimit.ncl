load "tests/getHistoryPath.ncl"

undef("saveSizeLimit")
procedure saveSizeLimit(historyRun:string,pattern:string)

;  inputs: historyRun, i.e. "coag04"
;          year, i.e. "1995"

local historyNames, rmin, cpi, rmassmin, readfile, writefile, hrc, nhrc, \
   allhist, time, nTimes, lev, nLevs, lat, nLats, lon, nLons, dustndbin, \
   dustrad, mixsulndbin, sulfrad, mixsul05a, mixsul10a, mixsul13a, mixsul30a, \
   dust05a, dust10a, dust13a, dust30a

begin

; Determine the path to history files
  historyPath = getHistoryPath(historyRun)
  
  historyPath = historyPath + "new/"
  
; Create array of files matching pattern
  allNames = systemfunc("ls " + historyPath + pattern)
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)
  
  nFiles=dimsizes(allNames)
  do ifile = 0, nFiles-1
    print("")
    print("file " + (ifile+1) + " of " + nFiles + " :")
    print("   " + (/ allNames(ifile) /))

    allhist = addfile(allNames(ifile), "w")
    
    ; Find the last time.
    time = allhist->time
    nTimes = dimsizes(time)
    lev = allhist->lev
    nLevs = dimsizes(lev)
    lat = allhist->lat
    nLats = dimsizes(lat)
    lon = allhist->lon
    nLons = dimsizes(lon)

    dustndbin=allhist->DUSTNDBIN ; ( time, lev, lat, lon, bins )
    dustrad=allhist->DUSTRAD ; ( bins )
;    dustdr=allhist->DUSTDR ; ( bins )
    mixsulndbin=allhist->MIXSULNDBIN ; ( time, lev, lat, lon, bins )
    sulfrad=allhist->SULFRAD  ; ( time, lev, lat, lon, bins )
    pursulndbin=allhist->PURSULNDBIN ; ( time, lev, lat, lon, bins )
    purerad=allhist->PURERAD  ; ( time, lev, lat, lon, bins )
    
    nBins=dimsizes(dustrad)

    ; get metadata from history
    mixsul05a = allhist->MIXSULND
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
;      print("bin "+iBin)

      if (.not.(isfilevar(allhist,"DUST05A"))) then
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
      end if
      
      if (.not.(isfilevar(allhist,"MIXSUL05A"))) then
        if (.not. mixThresh) then
          if (min(sulfrad(:, :, :, :, iBin)).ge.3.e-7) then
            mixThresh=True
          end if
        end if

        if (mixThresh) then
          mixsul05a(:, :, :, :) = mixsul05a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
          mixsul10a(:, :, :, :) = mixsul10a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
          mixsul13a(:, :, :, :) = mixsul13a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
          mixsul30a(:, :, :, :) = mixsul30a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)  
        else     
          mixsul05a(:, :, :, :) = mixsul05a(:, :, :, :) + where(sulfrad(:, :, :, :, iBin).ge.5e-8,   mixsulndbin(:, :, :, :, iBin), 0.0)
          mixsul10a(:, :, :, :) = mixsul10a(:, :, :, :) + where(sulfrad(:, :, :, :, iBin).ge.1e-7,   mixsulndbin(:, :, :, :, iBin), 0.0)
          mixsul13a(:, :, :, :) = mixsul13a(:, :, :, :) + where(sulfrad(:, :, :, :, iBin).ge.1.3e-7, mixsulndbin(:, :, :, :, iBin), 0.0)
          mixsul30a(:, :, :, :) = mixsul30a(:, :, :, :) + where(sulfrad(:, :, :, :, iBin).ge.3e-7,   mixsulndbin(:, :, :, :, iBin), 0.0)
        end if
      end if
      
      if (.not. pureThresh) then
        if (min(purerad(:, :, :, :, iBin)).ge.3.e-7) then
          pureThresh=True
        end if
      end if
      
      if (pureThresh) then
        pursul05a(:, :, :, :) = pursul05a(:, :, :, :) + pursulndbin(:, :, :, :, iBin)
        pursul10a(:, :, :, :) = pursul10a(:, :, :, :) + pursulndbin(:, :, :, :, iBin)
        pursul13a(:, :, :, :) = pursul13a(:, :, :, :) + pursulndbin(:, :, :, :, iBin)
        pursul30a(:, :, :, :) = pursul30a(:, :, :, :) + pursulndbin(:, :, :, :, iBin)  
      else
        pursul05a(:, :, :, :) = pursul05a(:, :, :, :) + where(purerad(:, :, :, :, iBin).ge.5e-8,   pursulndbin(:, :, :, :, iBin), 0.0)
        pursul10a(:, :, :, :) = pursul10a(:, :, :, :) + where(purerad(:, :, :, :, iBin).ge.1e-7,   pursulndbin(:, :, :, :, iBin), 0.0)
        pursul13a(:, :, :, :) = pursul13a(:, :, :, :) + where(purerad(:, :, :, :, iBin).ge.1.3e-7, pursulndbin(:, :, :, :, iBin), 0.0)
        pursul30a(:, :, :, :) = pursul30a(:, :, :, :) + where(purerad(:, :, :, :, iBin).ge.3e-7,   pursulndbin(:, :, :, :, iBin), 0.0)
      end if

      sec=sec2
      sec2=stringtointeger(systemfunc("date +%s"))
      sec=sec2-sec
      print("bin "+iBin+": "+sec+" sec")
    end do ; iBin
    
    delete(mixsulndbin)
    delete(sulfrad)
    
    delete(dustndbin)
    delete(dustrad)
    
    delete(pursulndbin)
    delete(purerad)
        
    delete(sec)
    delete(sec2)
    
    if (.not.(isfilevar(allhist,"MIXSUL05A"))) then
      print("writing output for MIXSUL...")
      allhist->MIXSUL05A=mixsul05a
      delete(mixsul05a)
      allhist->MIXSUL10A=mixsul10a
      delete(mixsul10a)
      allhist->MIXSUL13A=mixsul13a
      delete(mixsul13a)
      allhist->MIXSUL30A=mixsul30a
      delete(mixsul30a)
    end if
    
    if (.not.(isfilevar(allhist,"DUST05A"))) then
      print("writing output for DUST...")
      allhist->DUST05A=dust05a
      delete(dust05a)
      allhist->DUST10A=dust10a
      delete(dust10a)
      allhist->DUST13A=dust13a
       delete(dust13a)
     allhist->DUST30A=dust30a
      delete(dust30a)
    end if
    
    if (.not.(isfilevar(allhist,"PURSUL05A"))) then
      print("writing output for PURSUL...")
      allhist->PURSUL05A=pursul05a
      delete(pursul05a)
      allhist->PURSUL10A=pursul10a
      delete(pursul10a)
      allhist->PURSUL13A=pursul13a
      delete(pursul13a)
      allhist->PURSUL30A=pursul30a
      delete(pursul30a)
    end if
    
    print("cleaning up...")
    
    delete(allhist)

    delete(time)
    delete(nTimes)
    delete(lev)
    delete(nLevs)
    delete(lat)
    delete(nLats)
    delete(lon)
    delete(nLons)
    
    print("returning...")
   
  end do ; ifile
  
end

undef("fixDustND")
procedure fixDustND(historyRun:string,pattern:string)
begin

; Determine the path to history files
  historyPath = getHistoryPath(historyRun)
  
; Create array of files matching pattern
  allNames = systemfunc("ls " + historyPath + "all." + pattern)
  origNames = systemfunc("ls " + historyPath + "hist/" + pattern)
  
  nFiles=dimsizes(allNames)
  do ifile = 0, nFiles-1
    print("")
    print("file " + (ifile+1) + " of " + nFiles + " :")
    print("   " + (/ allNames(ifile) /))
    print("   " + (/ origNames(ifile) /))

    allhist = addfile(allNames(ifile), "w")
;    orighist = addfile(origNames(ifile), "r")
;    allhist->DUSTND=orighist->DUSTND
    delete(allhist->DUSTND@_FillValue)
  
  end do

end


; Main

;fixDustND("125expuT15","125expuT15.cam2.h3.1997-0[45]*")
;fixDustND("125expuT15","125expuT15.cam2.h3.1997-06-[01]*")

;saveSizeLimit("125expuT15","all.125expuT15.cam2.h3.1997-06-3*") 
;saveSizeLimit("125expuT15","all.125expuT15.cam2.h3.1997-0[78]*") 

;saveSizeLimit("125expuT15","all.125expuT15.cam2.h3.1997-0[45]*") 
;saveSizeLimit("125expuT15","all.125expuT15.cam2.h3.1997-06-[01]*") 

;saveSizeLimit("125smin","all.125smin.cam2.*") 
;saveSizeLimit("125expuT15b","all.125expuT15b.cam2.h3*") 
;saveSizeLimit("125expuT15b","all.125expuT15b.cam2.h3.1998-04-30*") 

;saveSizeLimit("nogeo400SO2","all.nogeo400SO2.cam2.*") 
saveSizeLimit("125smin400","all.125smin400.cam2.*") 
;saveSizeLimit("125sminPtbo2","all.125sminPtbo2.cam2.*") 
