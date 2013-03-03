load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("averageGrowingSeason")
procedure averageGrowingSeason(caseid:string)

begin
  path="$MODEL_PATH/ccsm/case/"
  filename1=path+"b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg/h0/b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg.cam2.h0.2013-07.nc"
  a1=addfile(filename1,"r")
  print("reading LANDFRAC")
  LANDFRAC=a1->LANDFRAC ; LANDFRAC(time, lat, lon)
  
  filename=path+caseid+"/proc/"+caseid+".cam.growseas.2013-2021.nc"
  print("adding "+filename)
  a=addfile(filename,"r")
  print("reading GROWSEAS")
  growseas=a->GROWSEAS ; TSMN(time, lat, lon)
  lat=a->lat
  lon=a->lon
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  
  fileout=path+caseid+"/proc/"+caseid+".cam.growseas.2014-2021avg.nc"
  print("rm "+fileout)
  b=systemfunc("rm "+fileout)
  print("creating "+fileout)
  c=addfile(fileout,"c")
    
  print("Averaging growseas over years 2014-2021")
  
  growseasavg=dim_avg_n_Wrap(growseas(1:8,:,:),0)
  print("Saving growseas to "+fileout)
  c->GROWSEAS=growseasavg 
  
  print("Done!")

  return
end

caseid="b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678"
averageGrowingSeason(caseid)

caseid="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg"
averageGrowingSeason(caseid)
