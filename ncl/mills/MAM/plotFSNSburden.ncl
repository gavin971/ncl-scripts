load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

runNamesConley = (/ "ccsm4_cam4_ocn" , "sad_trop_test" /)    
rootpath = "/glade/p/cesm/wawg/mmills/archive/"

do i = 0,1
  filepath=rootpath+runNamesConley(i)+"/atm/proc/h2tm/"
  filename=runNamesConley(i)+".cam.h2.tropicsavg.nc"
  fileid=addfile(filepath+filename,"r")
  date=fileid->date
  d1=ind(date.eq.19910501)
  d2=ind(date.eq.19930501)
  nTimes=d2-d1+1
  if (i.eq.0) then
    time=fileid->time(d1:d2)
    mass=new((/5,nTimes/),float)
    mass!1="time"
  end if
  print(i+" "+runNamesConley(i)+" d1:"+d1+" d2:"+d2+" nTimes:"+nTimes)
  print(date)
  mass(i,:)=fileid->VOLC_MASS_C(d1:d2)
  delete(date)
end do
          
runNamesMills = (/ "b.e11.B55TRW5CN.f19_g16.Pin16-18.dgnumhi1um.001", \
              "b.e11.B55TRW5CN.f19_g16.Pin20-30.dgnh1um.noqbo.001" , \
              "b.e11.B55TRW5CN.f19_g16.Pin17-25.dgnh1um.noqbo.001" /)
rootpath = "/glade/scratch/mmills/archive/"           
do i = 0,3
  filepath=rootpath+runNamesMills(i)+"/atm/proc/h3tm/"
  filename=runNamesMills(i)+".cam.h3.tropicsavg.nc"
  fileid=addfile(filepath+filename,"r")
  date=fileid->date
  d1=ind(date.eq.19910501)
  d2=ind(date.eq.19930501)
  nTimes=d2-d1+1
  print(i+" "+runNamesMills(i)+" d1:"+d1+" d2:"+d2+" nTimes:"+nTimes)
  mass(i+2,:)=fileid->TMso4_a1(d1:d2)
  mass(i+2,:)=mass(i+2,:)+fileid->TMso4_a2(d1:d2)
  delete(date)
end do

pdffile="MassBurden_Pinatubo"
wks = gsn_open_wks("pdf", "$IMAGE_PATH/"+pdffile)
res = True
res@gsnMaximize           = True
res@vpHeightF             = 0.5              ; change aspect ratio
res@vpWidthF              = 1.0
plot = gsn_csm_xy(wks, time, mass(:, :), res)
delete(wks)

