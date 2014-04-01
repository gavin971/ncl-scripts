load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "../nuclear/calculateMass.ncl"

procedure plotMass(basePath:string, historyPattern:string, runNames[*]:string, \
  tracerNames[*]:string, opt:logical)
  
; opt: optional argument, true if IC file is to be added via opt@icfile
  
local nRuns, nTracers, historyPath, history, nTimes, mass, wks, res, resP, i

begin

  nTracers = dimsizes(tracerNames)
    
  ; Calculate the tracer masses for all the runs.
  print("")
  print("")
  print("Run :  " + runNames)

  massfilepath = basePath + runNames + "/atm/proc/"
  massfile = runNames + "." + historyPattern + ".mass.nc"    
  print("mkdir -p "+massfilepath)
  a = systemfunc("mkdir -p "+massfilepath)
  historyPath = basePath + runNames + "/atm/hist/*." + historyPattern + ".*"
  print("calling calculateMass, historyPath="+historyPath)
  test = calculateMass(historyPath, tracerNames, opt)
  print("return from calculateMass")
  printVarSummary(test)
  time=test&time
  nt=dimsizes(time)


  print(" time@units:"+time@units+" nt="+nt)
  print(" dimsizes(test):"+dimsizes(test))
  ntmin=nt
  t1=0
  t2=ntmin
  t3=0

  print("t1="+t1+" t2="+t2+" t3="+t3+" ntmin="+ntmin)

  print("rm "+massfilepath+massfile)
  rm=systemfunc("rm "+massfilepath+massfile)
  print("creating "+massfilepath+massfile)
  mfileid = addfile(massfilepath+massfile, "c")
  ; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
  ;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
;      filedimdef(mfileid,"time",-1,True)

  so4_tot = new(nt,float)
  so4_tot(:)=0.0e0
  do j = 0, nTracers - 1
    k=j
    if (isStrSubset(tracerNames(j),"so4_a")) then
      so4_tot(:) = so4_tot(:) + (/test(j, :)/) ;- mass(i, j, 0)
    end if
    trcmass=test(j, :)
    trcmass!0="time"
    trcmass&time=time
    mfileid->$tracerNames(j)$=trcmass
  end do
  so4_tot!0="time"
  so4_tot&time=time
  mfileid->so4_tot=so4_tot
  delete(mfileid)

  delete(test)
  delete(trcmass)

  delete(so4_tot)
  delete(time)
  
end


rootpath="/glade/scratch/mmills/archive/"
;caseid = "b.e12.BRCP45W5CN.f19_g16.control.002"
caseid = "b.e12.BRCP45W5CN.f19_g16.geo10Tgpy2015TR.002"
                     
tracers=(/"SO2","H2SO4","so4_a2","so4_a1","so4_a3"/)
tracers@molwts=(/ 64.063799, 98.078479, 115.108999, 115.108999, 115.108999 /)
opt=False
opt@overwritemass=True
opt@lev1=1100.0
plotMass( rootpath, "h0", caseid, tracers, opt )
