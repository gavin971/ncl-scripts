
undef("saveSizeLimit")
procedure saveSizeLimit(historyRun:string,year:string)

;  inputs: historyRun, i.e. "coag04"
;          year, i.e. "1995"

local historyNames, rmin, cpi, rmassmin, readfile, writefile, hrc, nhrc, fc, sizehist, time, nTimes, lev, nLevs, lat, nLats, lon, nLons, dustndbin, dustrad, mixsulndbin, sulfrad, mixsul05a, mixsul10a, mixsul13a, mixsul30a, dust05a, dust10a, dust13a, dust30a, dustnd

begin
;  historyPattern=historyRun+"/dust-sulfur.cam2.h0."+year+"*.nc"
  sizePattern=historyRun+"/size.dust-sulfur.cam2.h0."+year+"*.nc"
; Create array of files matching pattern
;  historyNames = systemfunc("ls " + historyPattern)
  sizeNames = systemfunc("ls " + sizePattern)
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)
  
  nFiles=dimsizes(sizeNames)
  do ifile = 0, nFiles-1
    print("")
    print("    " + ifile + " :  " + (/ sizeNames(ifile) /))

    sizehist = addfile(sizeNames(ifile), "r")
    
    ; find corresponding history file
    fc=stringtochar(sizeNames(ifile))
    filen=historyRun+"/"+chartostring(fc(nhrc+5:))
    print(filen)
    history = addfile(filen,"r")
    
    filem=historyRun+"/limit."+chartostring(fc(nhrc+5:))
    print(filem)
    limit = addfile(filem,"c")
    
    ; Find the last time.
    time = sizehist->time
    nTimes = dimsizes(time)
    lev = sizehist->lev
    nLevs = dimsizes(lev)
    lat = sizehist->lat
    nLats = dimsizes(lat)
    lon = sizehist->lon
    nLons = dimsizes(lon)

    dustndbin=sizehist->DUSTNDBIN ; ( time, lev, lat, lon, bins )
    dustrad=sizehist->DUSTRAD ; ( bins )
;    dustdr=sizehist->DUSTDR ; ( bins )
    mixsulndbin=sizehist->MIXSULNDBIN ; ( time, lev, lat, lon, bins )
    sulfrad=sizehist->SULFRAD  ; ( time, lev, lat, lon, bins )
    
    nBins=dimsizes(dustrad)

    ; get metadata from history
    mixsul05a = history->MIXSULND
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
    
    do iBin = 0, nBins-1
      print("bin "+iBin)
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
    end do ; iBin
    
    limit->P0=history->P0
    limit->hyai=history->hyai
    limit->hyam=history->hyam 
    limit->hybm=history->hybm  
    limit->hybi=history->hybi 
    
    limit->MIXSUL05A=mixsul05a
    limit->MIXSUL10A=mixsul10a
    limit->MIXSUL13A=mixsul13a
    limit->MIXSUL30A=mixsul30a
    
    limit->DUST05A=dust05a
    limit->DUST10A=dust10a
    limit->DUST13A=dust13a
    limit->DUST30A=dust30a
    limit->DUSTND=dustnd
    
    delete(sizehist)
    delete(fc)
    delete(filen)
    delete(history)    
    delete(filem)
    delete(limit)

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
   
  end do    
end

; Main
;saveSizeLimit("coag04","1995-04-01")
;saveSizeLimit("coag04-2","1996-12-[23]")
saveSizeLimit("coag04-2","1996-12-31")
saveSizeLimit("coag04-2","1997")
