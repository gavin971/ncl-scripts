load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;caseid="b.e11.B1850W5CN.f19_g16.test_lim.001"
;filter=".h1.001"

rootdir="/glade/p/cesm/wawg/waccm/"
run_name = getenv("caseid")

if (ismissing(run_name)) then
  print("set environment variable caseid before running")
else
  filter = getenv("filter")
  if (ismissing(filter)) then
    filter=".h1."
    print("environment variable filter not found, defaulting to "+filter)
  end if
  outfilter = getenv("outfilter")
  if (ismissing(outfilter)) then
    outfilter = ""
  end if
  
  outfilename=rootdir+run_name+"/atm/proc/"+run_name+"."+outfilter+".Urms.nc"
  print("creating output file: "+outfilename)
  outfile=addfile(outfilename,"c")

  fils = systemfunc("ls "+rootdir+run_name+"/atm/hist/*"+filter+"*.nc")
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
