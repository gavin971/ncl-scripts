
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

rootpath = "/glade/scratch/mmills/archive/"

caseid = (/ "b.e12.B55TRW5CN.f19_g16.Pin20Tg.18-22km.1col.QBOon.001", \
         "b.e12.B55TRW5CN.f19_g16.Pin20Tg.18-22km.0-18N.96-118E.QBOon.001", \
         "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.1col.QBOon.001", \
         "b.e12.B55TRW5CN.f19_g16.Pin20Tg.20-27km.0-18N.96-118E.QBOon.001" /)
         
nRuns = dimsizes(caseid)

do i = 0,nRuns-1
  filepath=rootpath+caseid(i)+"/atm/proc/h0gm/"+caseid(i)+".cam.h0.globavg.nc"
  print(i+" "+filepath)
  a=addfile(filepath,"r")
  if (i.eq.0) then
    time=a->time
    nTimes=dimsizes(time)
    months=ispan(6,nTimes+5,1)
    months@long_name="Months since Jan 1, 1991"
    aodvis=new((/nRuns,nTimes/),"float")
    aodnir=new((/nRuns,nTimes/),"float")
    aoduv=new((/nRuns,nTimes/),"float")
    aodvisnir=new((/nRuns,nTimes/),"float")
    aoduvvis=new((/nRuns,nTimes/),"float")
    aoduvnir=new((/nRuns,nTimes/),"float")
  end if
  aodvis(i,:)=a->AODVISst
  aodnir(i,:)=a->AODNIRst
  aoduv(i,:)=a->AODUVst
  aodvisnir(i,:)=aodvis(i,:)/aodnir(i,:)
  aoduvvis(i,:)=aoduv(i,:)/aodvis(i,:)
  aoduvnir(i,:)=aoduv(i,:)/aodnir(i,:)
  delete(filepath)
  delete(a)
end do

aodvisnir@long_name="Vis/NIR"
aoduvvis@long_name="UV/Vis"
aoduvnir@long_name="UV/NIR"

pdfpath="$IMG_PATH/"
a = systemfunc("mkdir -p "+pdfpath)
pdffile=pdfpath+"ExtRatio"
print("Creating file "+pdffile+".pdf")
wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file

res                   = True
res@gsnPaperOrientation = "landscape"
res@gsnMaximize         = True        ; must include w/ Paper Orientation
res@vpHeightF             = 0.7              ; change aspect ratio
res@vpWidthF              = 1.0
res@xyLineThicknessF  = 4.0
res@xyMonoLineColor   = False
res@xyLineColors      = (/"black","red","blue","green"/)
res@xyMonoDashPattern = True
res@xyDashPattern = 0
res@

plot = gsn_csm_xy(wks, months, aodvis, res)
plot = gsn_csm_xy(wks, months, aodnir, res)
plot = gsn_csm_xy(wks, months, aoduv, res)
plot = gsn_csm_xy(wks, months, aodvisnir, res)
plot = gsn_csm_xy(wks, months, aoduvvis, res)
plot = gsn_csm_xy(wks, months, aoduvnir, res)
