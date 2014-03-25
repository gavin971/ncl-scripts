load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

rootpath="/glade/scratch/mmills/archive/"

caseid = (/ "b.e12.BRCP45W5CN.f19_g16.hifreq.001", \
            "b.e12.BRCP45W5CN.f19_g16.geo01Tgpy.001", \
            "b.e12.BRCP45W5CN.f19_g16.geo02Tgpy.001", \
            "b.e12.BRCP45W5CN.f19_g16.geo04Tgpy.001" /)

nCases = dimsizes(caseid)

do i = 0, nCases-1
  filepath=rootpath+caseid(i)+"/atm/proc/"+caseid(i)+"_RESTOM.nc"
  print(filepath)
  a=addfile(filepath,"r")
  r1=a->restom_ann_w
  nYears=dimsizes(r1)
  if (i.eq.0) then
    year=a->year
    print(year)
    restom=new((/nCases,nYears/),typeof(r1))
    restomc=restom
    flnt=restom
    fsnt=restom
    flntc=restom
    fsntc=restom
    restom_diff=new((/nCases-1,nYears/),typeof(r1))
    restomc_diff=restom_diff
  end if
  restom(i,0:nYears-1)=r1
  restomc(i,0:nYears-1)=a->restom_ann_w_c
  flntc(i,0:nYears-1)=a->flntc_ann_w
  fsntc(i,0:nYears-1)=a->fsntc_ann_w
  flnt(i,0:nYears-1)=a->flnt_ann_w
  fsnt(i,0:nYears-1)=a->fsnt_ann_w
  if (i.ne.0) then
    restom_diff(i-1,:)=restom(i,:)-restom(0,:)
    restomc_diff(i-1,:)=restomc(i,:)-restomc(0,:)
  end if
  delete(r1)
end do

print(restom)
print(restomc)
print(restom_diff)
print(restomc_diff)

image_path=getenv("IMAGE_PATH")
pdffile=image_path+"/RESTOM_diff"
print("Creating PDF file "+pdffile+".pdf")
wks   = gsn_open_wks ("pdf", pdffile)
  
plot = new(6, graphic)

resP                       = True
resP@gsnMaximize           = True            ; expand plot size to maximum
resP@gsnPaperOrientation   = "landscape"

res = True
res@gsnDraw               = False
res@gsnFrame           = False   ; Don't advance the frame

res@xyMonoDashPattern = True
res@xyMonoLineColor = False
res@xyLineColors = (/"green","blue","red"/)
res@xyLineThicknessF = 3

res@tiMainString = "RESTOM, expt-cntl"
plot(0) = gsn_csm_xy (wks,year,restom_diff,res)

res@tiMainString = "Clearsky RESTOM, expt-cntl"
plot(3) = gsn_csm_xy (wks,year,restomc_diff,res)

delete(res@xyLineColors)
res@xyLineColors = (/"black","green","blue","red"/)
res@tiMainString = "Clearsky LW net, TOM"
plot(4) = gsn_csm_xy (wks,year,flntc,res)

res@tiMainString = "Clearsky SW net, TOM"
plot(5) = gsn_csm_xy (wks,year,fsntc,res)

res@tiMainString = "Longwave net, Top Of Model"
plot(1) = gsn_csm_xy (wks,year,flnt,res)

res@tiMainString = "Shortwave net, TOM"
plot(2) = gsn_csm_xy (wks,year,fsnt,res)

gsn_panel(wks, plot, (/ 2, 3 /), resP)