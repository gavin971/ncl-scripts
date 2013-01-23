load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("surfIceAreaTseries")
procedure surfIceAreaTseries(caseid:string,path:string)

begin

  print(" ")
  print((/caseid/))

  threshold=0.9 ; fractional threshold for land, ice, and ocean masks

  npaths=dimsizes(path)
  
  historyNames = systemfunc("find " + path + " -name *.nc")
;  print((/historyNames/))
  nfiles=dimsizes(historyNames)

  lb=21
  latbounds=(/(/-90,90/),(/lb,90/),(/-90,-lb/)/)
  a=addfile(historyNames(0),"r")
  gw=a->gw
  lat=a->lat
  icefrac=a->ICEFRAC
  dsize=dimsizes(icefrac)
  nlat=dsize(1)
  nlon=dsize(2)

  wgt=new((/nlat,nlon/), typeof(gw))           
  wgt=conform(wgt,gw,0)
  wgt!0=gw!0
  wgt&lat=gw&lat 
    
  icefracave=new((/nfiles,3/), typeof(icefrac))


  do ifile=0,nfiles-1
    
    delete(a)
 ;   print((/historyNames(ifile)/))
    a=addfile(historyNames(ifile),"r")

    icefrac=a->ICEFRAC
    date=a->date

    do i = 0,2
;      print(" Latitude bounds: "+latbounds(i,0)+" to "+latbounds(i,1))

      icefracave(ifile,i)=wgt_areaave2(icefrac(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)

    end do
    
    print(date+" "+icefracave(ifile,0)+" "+icefracave(ifile,1)+" "+icefracave(ifile,2))
  end do
    
  print(" ")
  print((/caseid/))
  
  y1=156+5-96
  y2=y1+19

  do i = 0,2
    print(" ")
    print(" Latitude bounds: "+latbounds(i,0)+" to "+latbounds(i,1))
    print("  ice area, 96-295: "+avg(icefracave(:,i))+", stddev="+stddev(icefracave(:,i)))
    print("  ice area, 161-170: "+avg(icefracave(y1:y2,i))+", stddev="+stddev(icefracave(y1:y2,i)))
  end do
    
end

surfIceAreaTseries("b40.1850.track1.2deg.wcm.007","/data/waccm/b40.1850.track1.2deg.wcm.007/h0.annual")
