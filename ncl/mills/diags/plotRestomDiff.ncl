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

;caseid = (/ "b.e12.BRCP45W5CN.f19_g16.hifreq.001", \
;            "b.e12.BRCP45W5CN.f19_g16.geo01Tgpy.001", \
;            "b.e12.BRCP45W5CN.f19_g16.geo02Tgpy.001", \
;            "b.e12.BRCP45W5CN.f19_g16.geo04Tgpy.001", \
;            "b.e12.BRCP45W5CN.f19_g16.geo10Tgpy.001" /)

;caseid = (/ "f.e12.FSDW5.f19_f19.geo00.01", \
;            "f.e12.FSDW5.f19_f19.geo01.01", \
;            "f.e12.FSDW5.f19_f19.geo02.01", \
;            "f.e12.FSDW5.f19_f19.geo04.01", \
;            "f.e12.FSDW5.f19_f19.geo10.01" /)

caseid = (/ "p.e12.PSDW5.f19_f19.geo00s00.01", \
            "p.e12.PSDW5.f19_f19.geo00s01.01", \
            "p.e12.PSDW5.f19_f19.geo00s02.01", \
            "p.e12.PSDW5.f19_f19.geo00s04.01", \
            "p.e12.PSDW5.f19_f19.geo00s10.01" /)

nCases = dimsizes(caseid)

doTropFlux=False
print("nCases="+nCases)
do i = 0, nCases-1
  filepath=rootpath+caseid(i)+"/atm/proc/"+caseid(i)+"_RESTOM.nc"
  print(filepath)
  a=addfile(filepath,"r")
  r1=a->restom_ann_w
  nYears=dimsizes(r1)
  if (i.eq.0) then
    year=a->year
    print(year)
    nY=nYears
    restom=new((/nCases,nYears/),typeof(r1))
    restomc=restom
    flnt=restom
    fsnt=restom
    flntc=restom
    fsntc=restom
    if (isfilevar(a,"flnr_ann_w")) then
      flnr=restom
      fsnr=restom
      restrop=restom
      doTropFlux=True
    end if
    restom_diff=new((/nCases-1,nYears/),typeof(r1))
    restomc_diff=restom_diff
    flnt_diff=restom_diff
    fsnt_diff=restom_diff
    flntc_diff=restom_diff
    fsntc_diff=restom_diff
    if (doTropFlux) then
      flnr_diff=restom_diff
      fsnr_diff=restom_diff
      restrop_diff=restom_diff
    end if
  end if
  restom(i,0:nYears-1)=r1
  restomc(i,0:nYears-1)=a->restom_ann_w_c
  flntc(i,0:nYears-1)=a->flntc_ann_w
  fsntc(i,0:nYears-1)=a->fsntc_ann_w
  flnt(i,0:nYears-1)=a->flnt_ann_w
  fsnt(i,0:nYears-1)=a->fsnt_ann_w
  if (doTropFlux) then
    flnr(i,0:nYears-1)=a->flnr_ann_w
    fsnr(i,0:nYears-1)=a->fsnr_ann_w
    restrop(i,0:nYears-1)=fsnr(i,0:nYears-1)-flnr(i,0:nYears-1)
  end if
  if (i.ne.0) then
    restom_diff(i-1,:)=restom(i,:)-restom(0,:)
    restomc_diff(i-1,:)=restomc(i,:)-restomc(0,:)
    flnt_diff(i-1,:)=flnt(i,:)-flnt(0,:)
    fsnt_diff(i-1,:)=fsnt(i,:)-fsnt(0,:)
    flntc_diff(i-1,:)=flntc(i,:)-flntc(0,:)
    fsntc_diff(i-1,:)=fsntc(i,:)-fsntc(0,:)
    if (doTropFlux) then
      flnr_diff(i-1,:)=flnr(i,:)-flnr(0,:)
      fsnr_diff(i-1,:)=fsnr(i,:)-fsnr(0,:)
      restrop_diff(i-1,:)=restrop(i,:)-restrop(0,:)
    end if
  end if
  delete(r1)
end do

fstr="%7.3f"
line ="--------------------------------------------------------------------------------------"
line2="======================================================================================"
print((/line2/))
print("")
do i = 0,nCases-1
  if (i.eq.0) then
    print((/line2/))
    print("Control run")
    print((/line/))
    print((/caseid(i)/))
    print("YEAR   RESTOM   RESTOMC  RESTROP  FLNT     FSNT     FLNTC    FSNTC    FLNR     FSNR")
    print(year(:)+"  "+sprintf(fstr,restom(i,:))+"  "+sprintf(fstr,restomc(i,:))+"  "+\
       sprintf(fstr,restrop(i,:))+"  "+sprintf(fstr,flnt(i,:))+"  "+\
       sprintf(fstr,fsnt(i,:))+"  "+sprintf(fstr,flntc(i,:))+"  "+\
       sprintf(fstr,fsntc(i,:))+"  "+sprintf(fstr,flnr(i,:))+"  "+sprintf(fstr,fsnr(i,:)))
    print((/line2/))
    print("")
    print((/line2/))
    print("Experiments - control run:")
  else
    print((/line/))
    print((/caseid(i)/))
    print("YEAR  dRESTOM  dRESTOMC dRESTROP  dFLNT    dFSNT    dFLNTC   dFSNTC   dFLNR    dFSNR")
    print(year(:)+"  "+sprintf(fstr,restom_diff(i-1,:))+"  "+sprintf(fstr,restomc_diff(i-1,:))+"  "+\
       sprintf(fstr,restrop_diff(i-1,:))+"  "+sprintf(fstr,flnt_diff(i-1,:))+"  "+\
       sprintf(fstr,fsnt_diff(i-1,:))+"  "+sprintf(fstr,flntc_diff(i-1,:))+"  "+\
       sprintf(fstr,fsntc_diff(i-1,:))+"  "+sprintf(fstr,flnr_diff(i-1,:))+"  "+sprintf(fstr,fsnr_diff(i-1,:)))
  end if
end do
print((/line2/))

if (nYears.gt.1) then
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
  res@xyLineColors = (/"green","blue","red","brown"/)
  res@xyLineThicknessF = 3

  res@tiMainString = "RESTOM, expt-cntl"
  plot(0) = gsn_csm_xy (wks,year,restom_diff,res)

  txres = True
  txres@txFontHeightF = 0.03                     ; font smaller. default big
  txres@txFontColor="green"
  text1=gsn_add_text(wks,plot(0),"1 Tg/yr",2009,0.6,txres)
  txres@txFontColor="blue"
  text2=gsn_add_text(wks,plot(0),"2 Tg/yr",2009,0.3,txres)
  txres@txFontColor="red"
  text4=gsn_add_text(wks,plot(0),"4 Tg/yr",2009,0.0,txres)
  txres@txFontColor="brown"
  text4=gsn_add_text(wks,plot(0),"10 Tg/yr",2009,-0.3,txres)

  res@tiMainString = "Clearsky RESTOM, expt-cntl"
  plot(3) = gsn_csm_xy (wks,year,restomc_diff,res)

  delete(res@xyLineColors)
  res@xyLineColors = (/"black","green","blue","red","brown"/)
  res@tiMainString = "Clearsky LW net, TOM"
  plot(4) = gsn_csm_xy (wks,year,flntc,res)

  res@tiMainString = "Clearsky SW net, TOM"
  plot(5) = gsn_csm_xy (wks,year,fsntc,res)

  res@tiMainString = "Longwave net, Top Of Model"
  plot(1) = gsn_csm_xy (wks,year,flnt,res)
  txres@txFontColor="black"
  text0=gsn_add_text(wks,plot(1),"control",2009,236.1,txres)

  res@tiMainString = "Shortwave net, TOM"
  plot(2) = gsn_csm_xy (wks,year,fsnt,res)

  gsn_panel(wks, plot, (/ 2, 3 /), resP)

  delete(wks)

  print("evince "+pdffile+".pdf")
  evincer=systemfunc("evince " + pdffile + ".pdf &")
  print(evincer)
end if
