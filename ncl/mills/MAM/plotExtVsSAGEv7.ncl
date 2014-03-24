
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

rootpath="/glade/scratch/mmills/archive/"
runs=(/"f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.QBOon.001",\
       "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.QBOon.001",\
       "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.hetchem.002",\
       "f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.a3dgnl44.001", \
       "f.e12.FSDW5.f19_f19.Pin15Tg.19-27km.0N-20N.120E.a3dgnl44.001"/)


;       "f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.10S-20N.zavg.QBOon.001",\
;       "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.10S-20N.zavg.QBOon.001"/)
       
nRuns=dimsizes(runs)

datadir="/glade/p/work/mmills/inputdata/sulf/"

tropicsfile="CESM_1991_1995_sad_V2_sum_gw_SO2_TROPZ_tropicsavg.nc"
NETfile="CESM_1991_1995_sad_V2_sum_gw_SO2_TROPZ_NETavg.nc"
SETfile="CESM_1991_1995_sad_V2_sum_gw_SO2_TROPZ_SETavg.nc"

a=addfile(datadir+tropicsfile,"r")
date=a->date
d1=ind(date.eq."19910615")
print("d1="+d1)
nMonths=dimsizes(date(d1:))
AODstr_trop=new((/nRuns+1,nMonths/),"float")
AODstr_NET=AODstr_trop
AODstr_SET=AODstr_trop
AODstr_trop(0,:)=a->aod532strat(d1:)
delete(a)
delete(date)

a=addfile(datadir+NETfile,"r")
AODstr_NET(0,:)=a->aod532strat(d1:)
delete(a)

a=addfile(datadir+SETfile,"r")
AODstr_SET(0,:)=a->aod532strat(d1:)
delete(a)
delete(d1)

do i = 0,nRuns-1
  n=i+1
  print("run "+n)
  tropicsfile=rootpath+runs(i)+"/atm/proc/h0tm/"+runs(i)+".cam.h0.tropicsavg.nc"
  a=addfile(tropicsfile,"r")
  date=a->date
  d1=ind(date.eq."19910701")
  print("d1="+d1)
  nTimes=min((/dimsizes(date(d1:)),nMonths/))
  AODstr_trop(n,0:nTimes-1)=a->AODVISst(d1:d1+nTimes-1)
  delete(a)
  delete(date)
  
  NETfile=rootpath+runs(i)+"/atm/proc/h0gm/"+runs(i)+".cam.h0.NETavg.nc"
  a=addfile(NETfile,"r")
  AODstr_NET(n,0:nTimes-1)=a->AODVISst(d1:d1+nTimes-1)
  delete(a)
  
  SETfile=rootpath+runs(i)+"/atm/proc/h0gm/"+runs(i)+".cam.h0.SETavg.nc"
  a=addfile(SETfile,"r")
  AODstr_SET(n,0:nTimes-1)=a->AODVISst(d1:d1+nTimes-1)
  delete(a)
  delete(d1)
end do

month1=1991+6.5/12.
month2=month1+(nMonths-1.)/12.
time=fspan(month1,month2,nMonths)
;  print(time)

; Plot the data.
pdfpath="$IMG_PATH/"
a = systemfunc("mkdir -p "+pdfpath)
datestr = systemfunc("date +%y%m%d_%H%M")
pdffile=pdfpath+"ExtVsSAGEv7_c"+datestr
print("Creating file "+pdffile+".pdf")
wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file


res                   = True
res@gsnPaperOrientation = "landscape"
res@gsnMaximize         = True        ; must include w/ Paper Orientation
res@vpHeightF             = 0.7              ; change aspect ratio
res@vpWidthF              = 1.0
res@xyLineThicknessF  = 4.0
res@xyMonoDashPattern = True
res@xyMonoLineColor   = False
;res@xyLineColors       = (/"black","red","blue","blue","blue","green"/)
res@xyLineColors       = (/"black","red","blue","orange","green","brown"/)
res@trXMinF = 1991.5
res@trXMaxF = time(nMonths-1)
res@tmXBMinorPerMajor = 11
res@tmXTMinorPerMajor = 11
res@tiXAxisString = "Year"
res@tiYAxisString = "Stratospheric AOD, visible band"

res@tiMainString     = "Tropics Average"
plot = gsn_csm_xy(wks, time, AODstr_trop, res)

res@tiMainString     = "Northern Extratropics Average"
plot = gsn_csm_xy(wks, time, AODstr_NET, res)

res@tiMainString     = "Southern Extratropics Average"
plot = gsn_csm_xy(wks, time, AODstr_SET, res)

delete(wks)
