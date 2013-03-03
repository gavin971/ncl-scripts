load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("calcGrowingSeason")
procedure calcGrowingSeason(caseid:string)

begin
  path="$MODEL_PATH/ccsm/case/"
  filename1=path+"b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg/h0/b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg.cam2.h0.2013-07.nc"
  a1=addfile(filename1,"r")
  print("reading LANDFRAC")
  LANDFRAC=a1->LANDFRAC ; LANDFRAC(time, lat, lon)
  
  filename=path+caseid+"/proc/"+caseid+".cam2.h3.surfvars.2013-2022.nc"
  print("adding "+filename)
  a=addfile(filename,"r")
  print("reading TSMN")
  TSMN=a->TSMN ; TSMN(time, lat, lon)
  lat=a->lat
  lon=a->lon
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  
  fileout=path+caseid+"/proc/"+caseid+".cam.growseas.2013-2022.nc"
  print("rm "+fileout)
  b=systemfunc("rm "+fileout)
  print("creating "+fileout)
  c=addfile(fileout,"c")
 
; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
  filedimdef(c,"time",-1,True)

  years=ispan(2013,2021,1)
  date=years*10000+0701
  date!0="time"
  c->date=date
    
  growseas=new((/9,nlat,nlon/),integer)
  growseas@units="days"
  growseas@long_name="Length of growing season"
  growseas!0="time"
  time=(years-2013)*365+181
  time@units="days since 2013-01-01 00:00:00"
  time@calendar="noleap"
  time@long_name="time"
  time!0="time"
  growseas&time=time
  growseas!1="lat"
  growseas&lat=lat
  growseas!2="lon"
  growseas&lon=lon
  
  growseas(:,:,:)=growseas@_FillValue
  
  tcrit=273.15 ; Freezing point in Kelvins
  
  growday=where(TSMN(:,:,:).gt.tcrit, 1, 0) ;create 0/1 array
  
  do year = 0,8
;   year=0
    print("year "+years(year))
    do lt=0,nlat-1
      if (lat(lt).ge.0) then ; NH growing season
        t1=year*365
        t2=(year+1)*365-1      
      else ; SH growing season
        t1=year*365+183
        t2=(year+1)*365+183-1      
      end if
      do ln=0,nlon-1
         if (LANDFRAC(0,lt,ln).gt.0.001) then
           r0 = dim_numrun_n(growday(t1:t2,lt,ln), 0, 0) ; count number of warm "runs"
           ndays=365
           do while (r0(ndays-1).eq.0.and.ndays.ge.1)
             ndays=ndays-1
             if (ndays.eq.0) then
               break
             end if
           end do
           growseas(year,lt,ln)=ndays
         end if  
       end do    
    end do
  end do
  
  print("Saving growseas to "+fileout)
  c->GROWSEAS=growseas
  
  fileout2=path+caseid+"/proc/"+caseid+".cam.growseas.2014-2018avg.nc"
  print("rm "+fileout2)
  b=systemfunc("rm "+fileout2)
  print("creating "+fileout2)
  d=addfile(fileout2,"c")
 
  print("Averaging growseas over years 2014-2018")
  
  growseasavg=dim_avg_n_Wrap(growseas(1:5,:,:),0)
  print("Saving growseas to "+fileout2)
  d->GROWSEAS=growseasavg 
  
  print("Done!")

  return
end

caseid="b40.rcp4_5.2deg.wcm.carma.bc0tg.006"
calcGrowingSeason(caseid)

caseid="b40.rcp4_5.2deg.wcm.carma.bc0tg.007"
calcGrowingSeason(caseid)

caseid="b40.rcp4_5.2deg.wcm.carma.bc0tg.008"
calcGrowingSeason(caseid)

caseid="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003"
calcGrowingSeason(caseid)

caseid="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.006"
calcGrowingSeason(caseid)

caseid="b40.rcp4_5.2deg.wcm.carma.bc5tgIP.007"
calcGrowingSeason(caseid)

path="$MODEL_PATH/ccsm/case/"

cmd="ncea -O "+path+"b40.rcp4_5.2deg.wcm.carma.bc0tg.00[678]/proc/*.cam.growseas.2014-2018avg.nc "\
   +path+"b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678/proc/b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678.cam.growseas.2014-2018avg.nc"
print(cmd)
a=systemfunc(cmd)

cmd="ncea -O "+path+"b40.rcp4_5.2deg.wcm.carma.bc5tgIP.00[367]/proc/*.cam.growseas.2014-2018avg.nc "\
   +path+"b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg/proc/b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg.cam.growseas.2014-2018avg.nc"
print(cmd)
b=systemfunc(cmd)

cmd="ncea -O "+path+"b40.rcp4_5.2deg.wcm.carma.bc0tg.00[678]/proc/*.cam.growseas.2013-2022.nc "\
   +path+"b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678/proc/b40.rcp4_5.2deg.wcm.carma.bc0tg.avg678.cam.growseas.2013-2022.nc"
print(cmd)
c=systemfunc(cmd)

cmd="ncea -O "+path+"b40.rcp4_5.2deg.wcm.carma.bc5tgIP.00[367]/proc/*.cam.growseas.2013-2022.nc "\
   +path+"b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg/proc/b40.rcp4_5.2deg.wcm.carma.bc5tgIP.avg.cam.growseas.2013-2022.nc"
print(cmd)
d=systemfunc(cmd)
