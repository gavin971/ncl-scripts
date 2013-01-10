

undef("surfaceTable")
procedure surfaceTable(caseid:string,climopath:string)

begin

  print(" ")
  print((/caseid/))

  threshold=0.9 ; fractional threshold for land, ice, and ocean masks
  
  a=addfile(climopath+"/"+caseid+"_ANN_climo.nc","r")
  gw=a->gw
  landfrac=a->LANDFRAC
  ocnfrac=a->OCNFRAC
  icefrac=a->ICEFRAC
  TS=a->TS
  PSL=a->PSL
  PRECT=a->PRECC ; get metadata
  PRECT=PRECT+a->PRECL
  
  dsize=dimsizes(TS)
  nlat=dsize(1)
  nlon=dsize(2)
  
  wgt=new((/nlat,nlon/), typeof(gw))           
  wgt=conform(wgt,gw,0)
  wgt!0=gw!0
  wgt&lat=gw&lat 
  
  landmask=landfrac
  ocnmask=ocnfrac
  icemask=icefrac
  landmask=where(landfrac.ge.threshold,1.0,0.0)
  ocnmask=where(ocnfrac.ge.threshold,1.0,0.0)
  icemask=where(icefrac.ge.threshold,1.0,0.0)
     
  TSlandfrac=TS ; get metadata
  TSocnfrac=TS  
  TSicefrac=TS 
  
  TSlandfrac=TSlandfrac*landmask
  TSocnfrac=TSocnfrac*ocnmask  
  TSicefrac=TSicefrac*icemask 
  
  latbounds=(/(/-90,90/),(/20,90/),(/-20,20/),(/-90,-20/)/)
  
  do i = 0,3
    print(" Latitude bounds: "+latbounds(i,0)+" to "+latbounds(i,1))
    
    wgtbounds=wgt({latbounds(i,0):latbounds(i,1)},:)
    
    landfracave=wgt_areaave2(landmask(0,{latbounds(i,0):latbounds(i,1)},:),wgtbounds,0)
    ocnfracave=wgt_areaave2(ocnmask(0,{latbounds(i,0):latbounds(i,1)},:),wgtbounds,0)
    icefracave=wgt_areaave2(icemask(0,{latbounds(i,0):latbounds(i,1)},:),wgtbounds,0)
    
    TS_LAND=wgt_areaave2(TSlandfrac(0,{latbounds(i,0):latbounds(i,1)},:),wgtbounds,0)/landfracave  
    TS_OCN=wgt_areaave2(TSocnfrac(0,{latbounds(i,0):latbounds(i,1)},:),wgtbounds,0)/ocnfracave  
    
    if (icefracave.gt.0.0) then
     TS_ICE=wgt_areaave2(TSicefrac(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)/icefracave
    else
     TS_ICE=-999
    end if
    
    TS_AVE=wgt_areaave2(TS(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)
    PSL_AVE=wgt_areaave2(PSL(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)/100.
    PRECT_AVE=wgt_areaave2(PRECT(0,{latbounds(i,0):latbounds(i,1)},:),wgt({latbounds(i,0):latbounds(i,1)},:),0)
    PRECT_AVE=PRECT_AVE*1e3*60*60*24.
    
    SST=TS_OCN-273.15
    
    print("   TS_AVE ="+TS_AVE)
    print("   TS_LAND="+TS_LAND+", landfracave="+landfracave)
    print("   SST    ="+SST +", ocnfracave ="+ocnfracave)
    print("   PSL    ="+PSL_AVE)
    print("   PRECT  ="+PRECT_AVE)
    print("   TS_OCN ="+TS_OCN +", ocnfracave ="+ocnfracave)
    print("   TS_ICE ="+TS_ICE +", icefracave ="+icefracave)
    
    delete(wgtbounds)
  end do
  
end

surfaceTable("b40.1850.track1.2deg.wcm.007","/data/waccm/b40.1850.track1.2deg.wcm.007/climo/0096-0295")
surfaceTable("b40.1850.track1.1deg.006","/data/ccsm4/b40.1850.track1.1deg.006/climo/0800-0999")
surfaceTable("b40.1850.2deg.wset.002","/data/ccsm4/b40.1850.2deg.wset.002/climo/0020-0104")
