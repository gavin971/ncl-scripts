
maskFile="/glade/p/cesm/lmwg/diag/lnd_diag4.2/obs_data/clamp_data/" \
        + "observed/biomass/amazon_mask_1.9.nc"

maskFileID = addfile(maskFile,"r")
print("reading mask_amazon")
maskAmazon = maskFileID->mask_amazon

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
  cntlFire=exptid->COL_FIRE_CLOSS
  if (i.eq.0) then
    print("reading area")
    area=exptid->area
    areaAmazon=area*maskAmazon
    areaAmazonSum=sum(areaAmazon)
    print("Amazon area = "+areaAmazonSum)
    ds=dimsizes(exptFire)
    nYears=ds(0)
    exptFireAmazon = new((/nRuns,nYears/),"float")
    cntlFireAmazon = new((/nRuns,nYears/),"float")
  end if
  do j = 0,nYears-1
    print("area averaging year "+j)
    exptFireAmazon(i,j)=sum(exptFire(j,:,:)*areaAmazon)/areaAmazonSum
    cntlFireAmazon(i,j)=sum(cntlFire(j,:,:)*areaAmazon)/areaAmazonSum
  end do
  delete(exptid)
  delete(cntlid)
end do

do j = 0,nYears-1
  print(exptFireAmazon(0,j) + " " + exptFireAmazon(1,j) + " " + exptFireAmazon(2,j) + " " + \
        cntlFireAmazon(0,j) + " " + cntlFireAmazon(1,j) + " " + cntlFireAmazon(2,j))
end do
