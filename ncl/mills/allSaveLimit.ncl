
undef("saveSizeLimit")
procedure saveSizeLimit(historyRun:string,pattern:string)

;  inputs: historyRun, i.e. "coag04"
;          year, i.e. "1995"

local historyNames, rmin, cpi, rmassmin, readfile, writefile, hrc, nhrc, allhist, time, nTimes, lev, nLevs, lat, nLats, lon, nLons, dustndbin, dustrad, mixsulndbin, sulfrad, mixsul05a, mixsul10a, mixsul13a, mixsul30a, dust05a, dust10a, dust13a, dust30a, dustnd

begin

; Determine the path to history files
  if (isfilepresent("csm")) then
    historyPath="csm/"+historyRun+"/atm/hist/"
  else
    if (isfilepresent("waccm")) then
      historyPath="waccm/dust-sulfur/run/"+historyRun+"/"
    else
      if (isfilepresent("run")) then
        historyPath="run/"+historyRun+"/"
      else
        if (isfilepresent(historyRun)) then
          historyPath=historyRun+"/"
        end if
      end if
    end if
  end if
  
; Create array of files matching pattern
  allNames = systemfunc("ls " + historyPath + pattern)
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)
  
  nFiles=dimsizes(allNames)
  do ifile = 0, nFiles-1
    print("")
    print("file " + ifile + "of " + nFiles + " :  " + (/ allNames(ifile) /))

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
    dustnd = mixsul05a
            
    mixsul05a@long_name="mixed sulfates > 0.5 nm number density"
    mixsul10a@long_name="mixed sulfates > 1.0 nm number density"
    mixsul13a@long_name="mixed sulfates > 1.3 nm number density"
    mixsul30a@long_name="mixed sulfates > 3.0 nm number density"
    
    dust05a@long_name="dust > 0.5 nm number density"
    dust10a@long_name="dust > 1.0 nm number density"
    dust13a@long_name="dust > 1.3 nm number density"
    dust30a@long_name="dust > 3.0 nm number density"
    dustnd@long_name="dust number density"
    
    ; tally particle concentrations greater than .5, 1.0, 1.3. and 3 nm 
    ; NOTE: this tally could increase due to partial bin counts, not included yet...
    
    cross_thresh=False  ; Flag for when smallest sulfates are greater than largest threshold
    
    sec2=stringtointeger(systemfunc("date +%s"))
    print("looping bins 0 to "+nBins+"...")
    do iBin = 0, nBins-1
;      print("bin "+iBin)
      dustnd(:, :, :, :) = dustnd(:, :, :, :) + dustndbin(:, :, :, :, iBin)
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
      
      if (.not. cross_thresh) then
        if (min(sulfrad(:, :, :, :, iBin)).ge.3.e-7) then
          cross_thresh=True
        end if
      end if
      
      if (cross_thresh) then
        mixsul05a(:, :, :, :) = mixsul05a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
        mixsul10a(:, :, :, :) = mixsul10a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
        mixsul13a(:, :, :, :) = mixsul13a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)
        mixsul30a(:, :, :, :) = mixsul30a(:, :, :, :) + mixsulndbin(:, :, :, :, iBin)  
      else
        do iTime = 0, nTimes-1
          do iLev = 0, nLevs-1
            do iLat = 0, nLats-1
              do iLon = 0, nLons-1
                if (sulfrad(iTime, iLev, iLat, iLon, iBin).ge.5e-8) then
                  mixsul05a(iTime, iLev, iLat, iLon) = mixsul05a(iTime, iLev, iLat, iLon) + mixsulndbin(iTime, iLev, iLat, iLon, iBin)
                  if (sulfrad(iTime, iLev, iLat, iLon, iBin).ge.1e-7) then
                    mixsul10a(iTime, iLev, iLat, iLon) = mixsul10a(iTime, iLev, iLat, iLon) + mixsulndbin(iTime, iLev, iLat, iLon, iBin)
                    if (sulfrad(iTime, iLev, iLat, iLon, iBin).ge.1.3e-7) then
                      mixsul13a(iTime, iLev, iLat, iLon) = mixsul13a(iTime, iLev, iLat, iLon) + mixsulndbin(iTime, iLev, iLat, iLon, iBin)
                      if (sulfrad(iTime, iLev, iLat, iLon, iBin).ge.3.e-7) then
                        mixsul30a(iTime, iLev, iLat, iLon) = mixsul30a(iTime, iLev, iLat, iLon) + mixsulndbin(iTime, iLev, iLat, iLon, iBin)
                      end if
                    end if
                  end if  
                end if            
              end do ; iLon
            end do ; iLat
          end do ; iLev
        end do ; iTime
      end if
    sec=sec2
    sec2=stringtointeger(systemfunc("date +%s"))
    sec=sec2-sec
    print("bin "+iBin+": "+sec+" sec")
    end do ; iBin
    
    print("writing output for MIXSUL..."
    allhist->MIXSUL05A=mixsul05a
    allhist->MIXSUL10A=mixsul10a
    allhist->MIXSUL13A=mixsul13a
    allhist->MIXSUL30A=mixsul30a
    
    print("writing output for DUST..."
    allhist->DUST05A=dust05a
    allhist->DUST10A=dust10a
    allhist->DUST13A=dust13a
    allhist->DUST30A=dust30a
    allhist->DUSTND=dustnd
    
    delete(allhist)

    delete(time)
    delete(nTimes)
    delete(lev)
    delete(nLevs)
    delete(lat)
    delete(nLats)
    delete(lon)
    delete(nLons)
    
    delete(mixsul05a)
    delete(mixsul10a)
    delete(mixsul13a)
    delete(mixsul30a)
    delete(mixsulndbin)
    delete(sulfrad)
    
    delete(dust05a)
    delete(dust10a)
    delete(dust13a)
    delete(dust30a)
    delete(dustnd)
    delete(dustndbin)
    delete(dustrad)
    
    delete(sec)
    delete(sec2)
   
  end do    
end

; Main
;saveSizeLimit("nopur2_ncar","1995-07-[123]")
;saveSizeLimit("JPL06_ncar","1995-10-08")
;saveSizeLimit("JPL06_ncar","1995-10-2")
;saveSizeLimit("JPL06_ncar","1995-11")
;saveSizeLimit("expuT125E01","all.expuT125E01.cam2.h3.*.nc")

saveSizeLimit("expuT125E01","all.expuT125E01.cam2.h3.1997-10-17*") 
saveSizeLimit("expuT125E01","all.expuT125E01.cam2.h3.1997-10-2*") 
saveSizeLimit("expuT125E01","all.expuT125E01.cam2.h3.1997-11*") 
saveSizeLimit("expuT125E01","all.expuT125E01.cam2.h3.1997-12*") 
saveSizeLimit("expuT125E01","all.expuT125E01.cam2.h3.1998*") 
