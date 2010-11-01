load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("climoSootLoadCalc")
procedure climoSootLoadCalc(season:string, rootPath:string) 
begin 

  fileName=rootPath+"vgsoot08_"+season+"_climo.nc"
  a = addfile(fileName,"r")
  SOOTMD_col=a->SOOTMD_col ;(time, lat, lon) , units = "ng/m2"

  lat=a->lat
  lon=a->lon
  time=a->time
  gw=a->gw  ; gauss weights
  nlat=dimsizes(lat)
  
  do l1=0,nlat-1
    sootmass=dim_sum_Wrap(SOOTMD_col(time|0,lat|l1,lon|:))*gw(l1)
    
    sootmass=sootmass*510065600*1e6 ; Surface area of the Earth (m2)
    sootmass=sootmass/1e9 ; ng->g
    sootmass=sootmass/1e6 ; g->tonnes
    
    print(lat(l1)+" : "+sootmass)
    
  end do

end

season="ANN"
;season="DJF"
;season="JJA"

rootPath="/Volumes/Data/Models/waccm/run/carma-soot/vgsoot08/climo/2035-2064/"

climoSootLoadCalc(season,rootPath)
