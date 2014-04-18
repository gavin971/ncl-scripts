load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  
host=getenv("HOSTNAME")

if(isStrSubset(host,"pfe"))
  rootpath="/nobackupp8/mjmills2/archive/"
else
  rootpath="/glade/scratch/mmills/archive/"
end if

caseid = (/ "f.e12.FSDW5.f19_f19.geo00.02", \
            "f.e12.FSDW5.f19_f19.geo01.02", \
            "f.e12.FSDW5.f19_f19.geo02.02", \
            "f.e12.FSDW5.f19_f19.geo04.02", \
            "f.e12.FSDW5.f19_f19.geo10.02" /)

portid = (/ "p.e12.PSDW5.f19_f19.geo00s00.02", \
            "p.e12.PSDW5.f19_f19.geo00s01.02", \
            "p.e12.PSDW5.f19_f19.geo00s02.02", \
            "p.e12.PSDW5.f19_f19.geo00s04.02", \
            "p.e12.PSDW5.f19_f19.geo00s10.02" /)

burden =(/0., 1.,2.,4.,10./)

nCases = dimsizes(caseid)

doTropFlux=True
print("nCases="+nCases)
do i = 0, nCases-1
  filepath=rootpath+caseid(i)+"/atm/proc/"+caseid(i)+"_RESTOM.nc"
  portpath=rootpath+portid(i)+"/atm/proc/"+portid(i)+"_RESTOM.nc"
  i1=i+1
  print(i1+":"+filepath)
  a=addfile(filepath,"r")
  print(i1+":"+portpath)
  b=addfile(portpath,"r")
  r1=a->flnt_ann_w
  if (i.eq.0) then
    rf=new((/4,nCases/),typeof(r1))
    fln=rf
    fsn=rf
    
    rf_diff=rf
    fln_diff=rf
    fsn_diff=rf
  end if
  fln(0,i)=avg(r1)
  fsn(0,i)=avg(a->fsnt_ann_w)
  fln(1,i)=avg(a->flnr_ann_w)
  fsn(1,i)=avg(a->fsnr_ann_w)
  fln(2,i)=avg(b->flnt_ann_w)
  fsn(2,i)=avg(b->fsnt_ann_w)
  fln(3,i)=avg(b->flnr_ann_w)
  fsn(3,i)=avg(b->fsnr_ann_w)
  rf(:,i)=fsn(:,i)-fln(:,i)

  rf_diff(:,i)=rf(:,i)-rf(:,0)
  fln_diff(:,i)=fln(:,i)-fln(:,0)
  fsn_diff(:,i)=fsn(:,i)-fsn(:,0)
  delete(r1)
end do

print(burden+" "+fln(0,:)+" "+fsn(0,:)+" "+rf_diff(0,:))

image_path=getenv("IMAGE_PATH")
pdffile=image_path+"/GeoFluxes"
print("Creating PDF file "+pdffile+".pdf")
wks   = gsn_open_wks ("pdf", pdffile)

plot = new(3, graphic)

resP                       = True
resP@gsnMaximize           = True            ; expand plot size to maximum
resP@gsnPaperOrientation   = "landscape"

res = True
res@gsnDraw               = False
res@gsnFrame           = False   ; Don't advance the frame

res@xyMonoDashPattern = False
res@xyDashPatterns = (/ 0, 0, 1, 1 /)
res@xyMonoLineColor = False
res@xyLineColors = (/"blue","red","blue","red"/)
res@xyLineThicknessF = 3
res@trYMaxF = 0.5
res@trYMinF = -5.0
res@tiXAxisString="Injection rate (Tg SO2/year)"
res@tiYAxisString="Flux change (W/m2)"

res@tiMainString = "Net SW - Net LW, Geoeng-control"
plot(0) = gsn_csm_xy (wks,burden,rf_diff,res)

res@tiMainString = "Net Shortwave, Geoeng-control"
plot(1) = gsn_csm_xy (wks,burden,fsn_diff,res)

res@tiMainString = "Net Longwave, Geoeng-control"
plot(2) = gsn_csm_xy (wks,burden,fln_diff,res)

txres = True
txres@txFontHeightF = 0.037                     ; font smaller. default big
txres@txFontColor="blue"
text1=gsn_add_text(wks,plot(0),"Top of model",7.,-0.3,txres)
txres@txFontColor="red"
text2=gsn_add_text(wks,plot(0),"Tropopause",3.,-2.3,txres)
txres@txFontColor="black"
text3=gsn_add_text(wks,plot(0),"Solid: SD-WACCM",5.,-3.5,txres)
text4=gsn_add_text(wks,plot(0),"Dashed: Offline RT",5.,-4.2,txres)

gsn_panel(wks, plot, (/ 1, 3 /), resP)

delete(wks)

print("evince "+pdffile+".pdf")
evincer=systemfunc("evince " + pdffile + ".pdf &")
print(evincer)
