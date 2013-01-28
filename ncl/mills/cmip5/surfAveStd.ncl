load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("surfAveStd")
procedure surfAveStd(caseid:string,root:string,path:string)

begin

  print(" ")
  print((/caseid/))

  threshold=0.9 ; fractional threshold for land, ice, and ocean masks

  npaths=dimsizes(path)
  
  historyNames = systemfunc("find " + root + path + " -name *.nc")
  print((/historyNames/))
  nfiles=dimsizes(historyNames)

  lb=21
  latbounds=(/(/-90,90/),(/lb,90/),(/-lb,lb/),(/-90,-lb/)/)
  a=addfile(historyNames(0),"r")
  gw=a->gw
  lat=a->lat
  TS=a->TS
  dsize=dimsizes(TS)
  nlat=dsize(1)
  nlon=dsize(2)

  wgt=new((/nlat,nlon/), typeof(gw))           
  wgt=conform(wgt,gw,0)
  wgt!0=gw!0
  wgt&lat=gw&lat 
    
  TS_LAND=new((/nfiles,4/), typeof(TS))
  TS_OCN=TS_LAND 
  TS_AVE=TS_LAND
  PSL_AVE=TS_LAND
  PRECT_AVE=TS_LAND


  do ifile=0,nfiles-1
    
    delete(a)
    print((/historyNames(ifile)/))
    a=addfile(historyNames(ifile),"r")

    landfrac=a->LANDFRAC
    ocnfrac=a->OCNFRAC
    icefrac=a->ICEFRAC
    TS=a->TS
    PSL=a->PSL
    PRECT=a->PRECC ; get metadata
    PRECT=PRECT+a->PRECL

    landmask=landfrac
    ocnmask=ocnfrac
    icemask=icefrac
    landmask=where(landfrac.ge.threshold,1.0,0.0)
    ocnmask=where(ocnfrac.ge.threshold,1.0,0.0)
;    icemask=where(icefrac.ge.threshold,1.0,0.0)

    TSlandfrac=TS ; get metadata
    TSocnfrac=TS  
;    TSicefrac=TS 

    TSlandfrac=TSlandfrac*landmask
    TSocnfrac=TSocnfrac*ocnmask  
;    TSicefrac=TSicefrac*icemask 

;    TS_zm=dim_avg_Wrap(TS(0,:,:))

    do i = 0,3
;      print(" Latitude bounds: "+latbounds(i,0)+" to "+latbounds(i,1))

      landfracave=wgt_areaave2(landmask(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)
      ocnfracave=wgt_areaave2(ocnmask(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)
;      icefracave=wgt_areaave2(icemask(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)

      TS_LAND(ifile,i)=wgt_areaave2(TSlandfrac(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)/landfracave  
      TS_OCN(ifile,i)=wgt_areaave2(TSocnfrac(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)/ocnfracave  

;      if (icefracave.gt.0.0) then
;       TS_ICE=wgt_areaave2(TSicefrac(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)/icefracave
;      else
;       TS_ICE=-999
;      end if

      TS_AVE(ifile,i)=wgt_areaave2(TS(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)
      PSL_AVE(ifile,i)=wgt_areaave2(PSL(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)/100.
      PRECT_AVE(ifile,i)=wgt_areaave2(PRECT(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)

    end do
  end do
  
  delete(a)
  
  PRECT_AVE=PRECT_AVE*1e3*60*60*24.
  SST=TS_OCN-273.15
  print(" ")
  print((/caseid/))

  surfvars=new((/4,5,2/), typeof(TS_AVE))
;  TSgstd=TSgave
;  TSlave=TSgave
;  TSlstd=TSgave
;  SSTave=TSgave
;  SSTstd=TSgave
;  PSLave=TSgave
;  PSLstd=TSgave
;  PRECTave=TSgave
;  PRECTstd=TSgave
  
  do i = 0,3
    print(" ")
    print(" Latitude bounds: "+latbounds(i,0)+" to "+latbounds(i,1))
    surfvars(i,0,0)=avg(TS_AVE(:,i))
    surfvars(i,0,1)=stddev(TS_AVE(:,i))
    surfvars(i,1,0)=avg(TS_LAND(:,i))
    surfvars(i,1,1)=stddev(TS_LAND(:,i))
    surfvars(i,2,0)=avg(SST(:,i))
    surfvars(i,2,1)=stddev(SST(:,i))
    surfvars(i,3,0)=avg(PSL_AVE(:,i))
    surfvars(i,3,1)=stddev(PSL_AVE(:,i))
    surfvars(i,4,0)=avg(PRECT_AVE(:,i))
    surfvars(i,4,1)=stddev(PRECT_AVE(:,i))
    
    print("  TS_AVE: "+avg(TS_AVE(:,i))+", stddev="+stddev(TS_AVE(:,i)))
    print("  TS_LAND: "+avg(TS_LAND(:,i))+", stddev="+stddev(TS_LAND(:,i)))
    print("  SST: "+avg(SST(:,i))+", stddev="+stddev(SST(:,i)))
    print("  PSL_AVE: "+avg(PSL_AVE(:,i))+", stddev="+stddev(PSL_AVE(:,i)))
    print("  PRECT_AVE: "+avg(PRECT_AVE(:,i))+", stddev="+stddev(PRECT_AVE(:,i)))
  end do
  
  c = systemfunc("mkdir -p "+root+caseid+"/proc")
  delete(c)
  
  b=addfile(root+caseid+"/proc/"+caseid+".surfvars.nc","c")
  b->latbounds=latbounds
  b->surfvars=surfvars
  delete(b)    
end

;surfAveStd("b40.1850.track1.2deg.wcm.007","/data/waccm/","b40.1850.track1.2deg.wcm.007/h0.annual")

;surfAveStd("b40.1955-2005.2deg.wcm.avg","/data/waccm/","b40.1955-2005.2deg.wcm.00[234]/h0.annual")
;surfAveStd("b40.1850.track1.1deg.006","/data/ccsm4/","b40.1850.track1.1deg.006/h0.annual")
;surfAveStd("b40.1850.2deg.wset.002","/data/ccsm4/","b40.1850.2deg.wset.002/h0.annual")
surfAveStd("b40.1955-2005.2deg.wset.avg456","/data/ccsm4/","b40.1955-2005.2deg.wset.00[456]/h0.annual")

;surfAveStd("b40.20th.track1.1deg.avg567","/data/ccsm4/","b40.20th.track1.1deg.00[567]/h0.annual")
