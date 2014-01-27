load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; Some constants we may need:
pi=acos(0.0)*2.d0
Rearth  = 6.37122d6               ; Earth's radius (m)
Rearth  = Rearth * 1.d2          ; m -> cm
DG2RAD  = pi / 180.d0              ; rad deg^-1
Aearth  = 4.d0*pi*Rearth*Rearth ; Earth's surface area in cm2

print("Aearth="+Aearth+" cm2")

rootpath="/glade/p/work/mmills/jenglish/"
cases=(/ "geo-1-bhn", "geo-2-bhn", "geo-5-bhn", "geo-10-bhn" /)

nCases = dimsizes(cases)

;do i = 0,nCases-1
i=3

  filepath=rootpath+cases(i)+"/h0ann/"+cases(i)+".cam2.h0.avg.2015.nc"
  print("adding "+filepath)
  a=addfile(filepath,"w")
  z3=a->Z3 ; Geopotential Height in meters
  z3=z3*100.0 ; meters -> cm
  nlev=dimsizes(z3&lev)
  
  print("Calculating zdel...")
  zdel = z3 ; Geopotential Height in centimeters
  print("zdel: min="+min(zdel)+" max="+max(zdel))
  zdel(:,0:nlev-2,:,:)=z3(:,0:nlev-2,:,:)-z3(:,1:nlev-1,:,:)
  print("zdel: min="+min(zdel)+" max="+max(zdel))
  
  levcol=a->PURSULMD ; g/cm3
  levcol=levcol*zdel ; g/cm2/level
 
  print("Calculating column ...")
  column = dim_sum_Wrap(levcol(time|:, lat|:, lon|:, lev|:))
  column=column*32.064999/98.078479 ; g H2SO4 -> g S
  column@units="g S/cm2"
  column@long_name="PURSUL column mass density"
  column@method="column densities scaled by 32.064999/98.078479 to get g of S"
  column@program="calcGeoSBurdens.ncl by Mike Mills, 2014.01.27"
  print("Writing column ...")
  a->PURSULMDcol=column ; (time, lat, lon)
  
  gw=a->gw      ; latitude weights
  colave=wgt_areaave_Wrap(column,gw,1.0,1)
  colsum=colave*Aearth/1.e12 ; g->Tg
  print("colave="+colave+"g S/cm2, colsum="+colsum+" Tg S")
  
  print("Calculating grid ...")
  nLons=dimsizes(z3&lon)
  nLats=dimsizes(z3&lat)
  gwtot=sum(gw) ; latitude weights sum to 2
  
  AearthPerLon=Aearth/nLons ; cm2 per longitude wedge
  AreaPerGrid = AearthPerLon*gw/gwtot
  
  column2=column ; get metadata
  do j = 0, nLats-1
;    print(j+" "+AreaPerGrid(j))
    column2(:,j,:)=column2(:,j,:)*doubletofloat(AreaPerGrid(j))
  end do
  column2=column2/1e12
  column2@units="Tg S/grid box"
  column2@long_name="PURSUL column mass per grid box"
  
  a->PURSULMDcolgrd=column2
  
  delete(a)
  delete(zdel)
  delete(levcol)
  delete(column)    
  delete(z3)
  delete(nlev)
  
;end do
