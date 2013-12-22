load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
IMAGE_PATH=getenv("IMG_PATH")

rootpath="/glade/p/acd/mmills/run/"

Cases = (/ "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003_2013-2038", \
           "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.006_2013-2038", \
           "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.007_2013-2038", \)
           "b40.rcp4_5.2deg.wcm.carma.bc0tg.006_2013-2038", \
           "b40.rcp4_5.2deg.wcm.carma.bc0tg.007_2013-2038", \
           "b40.rcp4_5.2deg.wcm.carma.bc0tg.008_2013-2038" /)
               
nRuns = dimsizes(Cases)

do i = 0,nRuns-1
  exptFile=rootpath+Cases(i)+"/"+Cases(i)+"_ANN_ALL.nc"
  exptid = addfile(exptFile,"r")
  print(i+" reading COL_FIRE_CLOSS from "+Cases(i))
  exptFire=exptid->COL_FIRE_CLOSS
  if (i.eq.0) then
    print("reading area")
    area=exptid->area
    areaAmazonSum=sum(area({0:-10},{290:310}))
    print("Amazon area = "+areaAmazonSum)
    ds=dimsizes(exptFire)
    nYears=ds(0)
    exptFireAmazon = new((/nRuns,nYears/),"float")
    exptFireAmazon@long_name=exptFire@long_name
  end if
  do j = 0,nYears-1
;    print("area averaging year "+j)
    exptFireAmazon(i,j)=sum(exptFire(j,{0:-10},{290:310})*area({0:-10},{290:310}))
  end do
  delete(exptid)
end do

exptFireAmazon=exptFireAmazon*3600.*24./areaAmazonSum
exptFireAmazon@units="gC/m^2/d"

exptFireAmazonAvg = avg(exptFireAmazon(0:2,0:7))
cntlFireAmazonAvg = avg(exptFireAmazon(3:5,0:7))

print("experiment average:"+exptFireAmazonAvg)
print("   control average:"+cntlFireAmazonAvg)
print(exptFireAmazonAvg/cntlFireAmazonAvg)

pdffile="LandFireAmazon"
wks=gsn_open_wks ("pdf",IMAGE_PATH+"/nuclear/"+pdffile)

res=True
res@xyMonoDashPattern = True
res@xyDashPattern = 0
res@xyMonoLineColor    = False   ; different colors for different lines
res@xyLineColors       = (/"red","red","red","blue","blue","blue"/)
res@xyLineThicknessF = 3.0
res@tiXAxisString="Year"
res@tiYAxisString="gC/m^2/d"
res@tiMainString=exptFireAmazon@long_name
;res@trXMinF=1
res@trXMaxF=nYears

years=ispan(1,nYears,1)
printVarSummary(years)
printVarSummary(exptFireAmazon)

baseplot = gsn_csm_xy (wks, years, exptFireAmazon, res)
delete(wks)


