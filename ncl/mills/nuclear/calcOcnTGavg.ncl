caseid = "b40.rcp4_5.2deg.wcm.carma.bc5tgIP.006"

rootpath = "/glade/p/acd/mmills/run/"

filename = rootpath+caseid+"/ocn/proc/"+caseid+".pop.h.TEMP.2013-2038.nc"

print("reading "+filename)
fileid = addfile(filename,"r")

print("reading TEMP")
TEMP=fileid->TEMP ; [time | 312] x [z_t | 60] x [nlat | 384] x [nlon | 320]
printVarSummary(TEMP)
ds = dimsizes(TEMP)
nTimes = ds(0)
nLev = ds(1)
nLat = ds(2)
nLon = ds(3)

print("reading REGION_MASK")
REGION_MASK=fileid->REGION_MASK ; [nlat | 384] x [nlon | 320]
printVarSummary(REGION_MASK)

print("reading TAREA")
TAREA=fileid->TAREA ; [nlat | 384] x [nlon | 320]
printVarSummary(TAREA)

print("creating mask") ; remove land and marginal sea regions 
ocnmask = where(REGION_MASK > 0, 1.0, 0.0)
printVarSummary(ocnmask)

print("creating ocnarea") 
ocnarea = TAREA * ocnmask
printVarSummary(ocnarea)
ocnareasum = sum(ocnarea)
print("ocnareasum = "+ocnareasum)

TempGavg = new((/nTimes,nLev/), "double")
TempGavg!0=TEMP!0
TempGavg!1=TEMP!1
TempGavg&time=TEMP&time
TempGavg&z_t=TEMP&z_t
printVarSummary(TempGavg)

do t = 0,nTimes-1
  t1 = t + 1
  print(t1+" of "+nTimes)
  do k = 0,nLev-1
    Tocnarea = TEMP(t,k,:,:)*ocnarea(:,:)
    TempGavg(t,k) = sum(Tocnarea)
  end do
end do  
TempGavg = TempGavg/ocnareasum

fileout=rootpath+caseid+"/ocn/proc/"+caseid+".pop.h.TEMP.Gavg.2013-2038.nc"
del=systemfunc("rm "+fileout)
print("creating "+fileout)
fileoutid = addfile(fileout,"c")
print("writing coordinate variables")
fileoutid->time=fileid->time
fileoutid->moc_components=fileid->moc_components
fileoutid->z_t=fileid->z_t
fileoutid->z_t_150m=fileid->z_t_150m
fileoutid->z_w=fileid->z_w
fileoutid->z_w_top=fileid->z_w_top
fileoutid->moc_z=fileid->moc_z
fileoutid->dzw=fileid->dzw
print("writing TEMP")
fileoutid->TEMP=TempGavg
delete(fileid)
delete(fileoutid)

