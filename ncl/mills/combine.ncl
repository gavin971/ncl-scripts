undef("combine")
procedure combine(historyRun:string,allDir:string,year:string)

;  inputs: historyRun, i.e. "coag04"
;          allDir, i.e. "coag04all"
;          year, i.e. "1995"

local historyNames,rmin,cpi,rmassmin,readfile,writefile,hrc,nhrc,fc

begin
  historyPattern=historyRun+"/dust-sulfur.cam2.h0."+year+"*.nc"
; Create array of files matching pattern
  historyNames = systemfunc("ls " + historyPattern)
  hrc=stringtochar(historyRun)
  nhrc=dimsizes(hrc)
  
  nFiles=dimsizes(historyNames)
  do ifile = 0, nFiles-1
    print("")
    print("    " + ifile + " :  " + (/ historyNames(ifile) /))

    history = addfile(historyNames(ifile), "r")
    fc=stringtochar(historyNames(ifile))
    sizefilename=historyRun+"/size."+chartostring(fc(nhrc:))
    print(sizefilename)
    limitfilename=historyRun+"/limit."+chartostring(fc(nhrc:))
    print(limitfilename)
    sizefile=addfile(sizefilename,"r")   
    limitfile=addfile(limitfilename,"r")
    
    ; copy history file:
    allfilename=allDir+"/all."+chartostring(fc(nhrc:))
    print(allfilename)
    out = systemfunc("cp " + historyNames(ifile) + " " + allfilename )
    writefile=addfile(allfilename,"w")

;    writefile->bins         = sizefile->bins
    writefile->DM           = sizefile->DM
    writefile->DUST05A      = limitfile->DUST05A
    writefile->DUST10A      = limitfile->DUST10A
    writefile->DUST13A      = limitfile->DUST13A
    writefile->DUST30A      = limitfile->DUST30A
    writefile->DUSTDR       = sizefile->DUSTDR
    writefile->DUSTNDBIN    = sizefile->DUSTNDBIN
    writefile->DUSTRAD      = sizefile->DUSTRAD
    writefile->DUSTRADLOW   = sizefile->DUSTRADLOW
    writefile->DUSTRADUP    = sizefile->DUSTRADUP
    writefile->DUSTVOL      = sizefile->DUSTVOL
    writefile->MIXSUL05A    = limitfile->MIXSUL05A
    writefile->MIXSUL10A    = limitfile->MIXSUL10A
    writefile->MIXSUL13A    = limitfile->MIXSUL13A
    writefile->MIXSUL30A    = limitfile->MIXSUL30A
    writefile->MIXSULNDBIN  = sizefile->MIXSULNDBIN
    writefile->MIXSULVMR    = sizefile->MIXSULVMR
    writefile->PURSULVMR    = sizefile->PURSULVMR
    writefile->PRESSURE     = sizefile->PRESSURE
    writefile->RHOAIR       = sizefile->RHOAIR
    writefile->RMASS        = sizefile->RMASS
    writefile->RMASSUP      = sizefile->RMASSUP
    writefile->SOY          = sizefile->SOY
    if (isfilevar(history,"SULDEN")) then
      writefile->SULDEN=history->SULDEN    
    else
      writefile->SULDEN       = sizefile->SULDEN
    end if
    writefile->SULFRAD      = sizefile->SULFRAD
    
  end do
end

combine("coag04","coag04all","1995-04")
;combine("coag04-2","coag04all","1995-09")
