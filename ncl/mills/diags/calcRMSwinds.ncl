load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

caseid = getenv("caseid") ; e.g. b.e11.B1850W5CN.f19_g16.test_lim.001

if (ismissing(caseid)) then

  print("set environment variable caseid before running")
  
else
  
  rootdir = getenv("rootdir")
  if (ismissing(filter)) then
    rootdir="/glade/p/cesm/wawg/waccm/"
    print("environment variable rootdir not found, defaulting to "+filter)
  end if

  filter = getenv("filter") ; e.g. ".h1.001[5-9]"
  if (ismissing(filter)) then
    filter=".h1."
    print("environment variable filter not found, defaulting to "+filter)
  end if
  
  outfilter = getenv("outfilter") ; e.g. "h1.0015-0019"
  if (ismissing(outfilter)) then
    outfilter = ""
    print("environment variable outfilter not found")
  end if
  
  outfilepath=rootdir+caseid+"/atm/proc/"
  a=systemfunc("mkdir -p "+outfilepath)
  outfilename=outfilepath+caseid+"."+outfilter+".Urms.nc"
  print("creating output file: "+outfilename)
  outfile=addfile(outfilename,"c")

  fils = systemfunc("ls "+rootdir+caseid+"/atm/hist/*"+filter+"*.nc")
  f1   = addfile(fils(0),"r")
  outfile->hyam=f1->hyam
  outfile->hybm=f1->hybm
  outfile->hyai=f1->hyai
  outfile->hybi=f1->hybi
  outfile->ilev=f1->ilev
  outfile->P0=f1->P0
  outfile->slat=f1->slat
  outfile->slon=f1->slon
  outfile->gw=f1->gw
  
  print("adding files:")
  print(fils)
  f    = addfiles(fils, "r")
  U    = addfiles_GetVar(f,fils,"U")   ; float U(time, lev, lat, lon) 
  printVarSummary(U)
  
  time=U&time
  lon =U&lon
;  dsizes=dimsizes(U)
  nTimes=(/dimsizes(time)/)
  nLons=(/dimsizes(lon)/)
  print("nTimes="+nTimes)
  print("nLons ="+nLons)
  
  print("calculating zonal mean...")
  U_zm = dim_avg_Wrap(U)
  
  Ubar = U
  
  print("calculating Ubar...")
  do iLon=0,nLons-1
    print("iLon="+iLon)
    Ubar(:,:,:,iLon)=U(:,:,:,iLon)-U_zm(:,:,:)
  end do
  print("squaring Ubar")
  Ubar=Ubar*Ubar
  
  print("calculating temporal mean...")
  Urms=dim_avg_Wrap(Ubar(lev|:, lat|:, lon|:, time|:))
  
  print("calculating square root...")
  Urms=sqrt(Urms)
  
  print("saving Urms...")
  outfile->Urms=Urms
  
  delete(outfile)
end if
