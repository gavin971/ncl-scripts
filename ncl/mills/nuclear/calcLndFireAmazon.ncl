
maskFile="/glade/p/cesm/lmwg/diag/lnd_diag4.2/obs_data/clamp_data/" \
        + "observed/biomass/amazon_mask_1.9.nc"

;maskFileID = addfile(maskFile,"r")
;print("reading mask_amazon")
;maskAmazon = maskFileID->mask_amazon

;maskAmazon=where(maskAmazon.lt.1.0, 0.0, maskAmazon) 

rootpath="/glade/p/acd/mmills/run/"

exptCases = (/ "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.003_2013-2038", \
               "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.006_2013-2038", \
               "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.007_2013-2038" /)

cntlCases = (/ "b40.rcp4_5.2deg.wcm.carma.bc0tg.006_2013-2038", \
               "b40.rcp4_5.2deg.wcm.carma.bc0tg.007_2013-2038", \
               "b40.rcp4_5.2deg.wcm.carma.bc0tg.008_2013-2038" /)
               
nRuns = dimsizes(exptCases)

do i = 0,nRuns-1
  exptFile=rootpath+exptCases(i)+"/"+exptCases(i)+"_ANN_ALL.nc"
  cntlFile=rootpath+cntlCases(i)+"/"+cntlCases(i)+"_ANN_ALL.nc"
  exptid = addfile(exptFile,"r")
  cntlid = addfile(cntlFile,"r")
  print(i+" reading COL_FIRE_CLOSS from "+exptCases(i))
  exptFire=exptid->COL_FIRE_CLOSS
  print(i+" reading COL_FIRE_CLOSS from "+cntlCases(i))
  cntlFire=cntlid->COL_FIRE_CLOSS
  if (i.eq.0) then
    print("reading area")
    area=exptid->area
;    areaAmazon=area*maskAmazon
    areaAmazonSum=sum(area({0:-10},{290:310}))
    print("Amazon area = "+areaAmazonSum)
    ds=dimsizes(exptFire)
    nYears=ds(0)
    exptFireAmazon = new((/nRuns,nYears/),"float")
    cntlFireAmazon = new((/nRuns,nYears/),"float")
  end if
  do j = 0,nYears-1
;    print("area averaging year "+j)
    exptFireAmazon(i,j)=sum(exptFire(j,{0:-10},{290:310})*area({0:-10},{290:310}))
    cntlFireAmazon(i,j)=sum(cntlFire(j,{0:-10},{290:310})*area({0:-10},{290:310}))
  end do
  delete(exptid)
  delete(cntlid)
end do

exptFireAmazon=exptFireAmazon*3600.*24./areaAmazonSum
cntlFireAmazon=cntlFireAmazon*3600.*24./areaAmazonSum

do j = 0,nYears-1
  print(exptFireAmazon(0,j) + " " + exptFireAmazon(1,j) + " " + exptFireAmazon(2,j) + " " + \
        cntlFireAmazon(0,j) + " " + cntlFireAmazon(1,j) + " " + cntlFireAmazon(2,j))
end do

exptFireAmazonAvg  = avg(exptFireAmazon(0,0:7))
exptFireAmazonAvg1 = avg(exptFireAmazon(1,0:7))
exptFireAmazonAvg2 = avg(exptFireAmazon(2,0:7))
cntlFireAmazonAvg  = avg(cntlFireAmazon(:,0:7))
cntlFireAmazonAvg0  = avg(cntlFireAmazon(0,0:7))
cntlFireAmazonAvg1  = avg(cntlFireAmazon(1,0:7))
cntlFireAmazonAvg2  = avg(cntlFireAmazon(2,0:7))

print("experiment run 1:"+exptFireAmazonAvg)
print("experiment run 2:"+exptFireAmazonAvg1)
print("experiment run 3:"+exptFireAmazonAvg2)
print(" control average:"+cntlFireAmazonAvg)
print("   control run 1:"+cntlFireAmazonAvg0)
print("   control run 2:"+cntlFireAmazonAvg1)
print("   control run 3:"+cntlFireAmazonAvg2)
print(exptFireAmazonAvg/cntlFireAmazonAvg)
print(exptFireAmazonAvg1/cntlFireAmazonAvg)
print(exptFireAmazonAvg2/cntlFireAmazonAvg)
