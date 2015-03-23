load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  res=True
  
  caseid = (/  \
             "b.e13.BHISTW5TCN.f19_g16.beta17_stratmam.001", \
             "b.e13.BHISTW5TCN.f19_g16.beta17_stratmam.001", \
             "b.e13.BHISTW5TCN.f19_g16.beta17_stratmam.001" \
               /)  
  yearString = (/ "1992-1993", "1996-1997", "1989-1990" /)
  res@xyMonoLineColor = False
  res@xyLineColors    = (/"red","green3", "black"/)
  res@xyMonoDashPattern = False
  res@xyDashPatterns = (/0,0,0,0/)
  
  nCases = dimsizes(caseid)
  
  do i = 0, nCases-1

    rootpath = "/glade/scratch/mmills/archive/"
    exists =isfilepresent(rootpath+caseid(i)+"/atm/proc/SizeDistN2O_"+yearString(i)+"_NH.nc")
    if (.not.exists) then
      print("Missing: "+rootpath+caseid(i)+"/atm/proc/h0zm/SizeDistN2O_"+yearString(i)+"_NH.nc")
      rootpath="/glade/p/cesm/wawg_dev/mmills/archive/"
    end if
    print("Rootpath: "+rootpath+caseid+"/atm/proc/h0zm/")
  
    a=addfile(rootpath+caseid(i)+"/atm/proc/SizeDistN2O_"+yearString(i)+"_NH.nc","r")
    if (i.eq.0) then
      DpMean=a->diameter
      nBins=dimsizes(DpMean)
      nB=nBins-1
      N2Obounds=a->N2Oppb
      nRanges=dimsizes(N2Obounds)-1
      dVdlog10Dp_avg=new((/nCases,nRanges,nBins/),"float")
      dVdlog10Dp_avg!2="diameter"
      dVdlog10Dp_avg&diameter=DpMean
      dVdlog10Dp_std=dVdlog10Dp_avg
    else  
      nB=dimsizes(a->diameter)
      nB=nB-1
    end if
    
    dVdlog10Dp_avg(i,:,:nB)=a->dVdlog10Dp_avg
    dVdlog10Dp_std(i,:,:nB)=a->dVdlog10Dp_std
  
  end do
  
  datestr = systemfunc("date +%y%m%d_%H%M")
              
  pdffile="SizeDistN2O_multi_"+caseid(0)+"_"+yearString(0)+"_NH_c"+datestr
  wks   = gsn_open_wks ("pdf", "$IMAGE_PATH/"+pdffile)                   ; open workstation
  res@vpHeightF             = 0.4              ; change aspect ratio
  res@vpWidthF              = 0.7
  res@gsnFrame               = False                     ; don't draw yet
  res@gsnDraw                = False                     ; don't advance frame
  res@xyXStyle="Log"
  res@trXMinF=0.05
  res@trXMaxF=5.0
  res@trYMinF=-0.2
  res@trYMaxF=3.1
  res@xyLineThicknessF = 4.0
  res@tiXAxisString = "Particle Diameter, Dp, ~F33~m~FR~m"
  res@tiYAxisString = "dV/dLog~B~10~N~Dp, normalized"
  res@tmXBMinorPerMajor = 8
  res@tmXTMinorPerMajor = 8
  res@tmXBMajorOutwardLengthF= -0.021
  res@tmXBMajorLengthF= -0.021
  res@tmXBMinorOutwardLengthF= -0.01
  res@tmXBMinorLengthF= -0.01
  res@tmYLMajorLengthF= 0.021
  res@tmYLMajorOutwardLengthF= 0.021
  res@tmYLMinorPerMajor = 3
  res@tmYRMinorPerMajor = 3
  res@tmYLMinorOn=True
  res@tmYLMinorLengthF= 0.008
  res@tmYLMinorOutwardLengthF= 0.008
  
  polyres                   = True                       ; marker resources
  polyres@gsMarkerIndex     = 1                          ; index for circle
  polyres@gsMarkerSizeF     = .02                        ; size 
  polyres@gsLineThicknessF  = 2
  
  dfac =1.05
  
  do k = 1, nRanges
    j=nRanges-k
    res@tiMainString = N2Obounds(j+1)+" < N2O < "+N2Obounds(j)
  
    plot = gsn_csm_xy (wks, DpMean, dVdlog10Dp_avg(:,j,:), res)

  ;************************************************
  ; add error bars
  ;************************************************ 
    error_bar = new(nCases*nBins,graphic)
    error_bar_top = new(nCases*nBins,graphic)
    error_bar_bot = new(nCases*nBins,graphic)
    do c = 0, nCases-1 
      polyres@gsLineColor = res@xyLineColors(c)

      do i = 0, nBins-1
        top=dVdlog10Dp_avg(c,j,i)+dVdlog10Dp_std(c,j,i)
        bottom=max((/dVdlog10Dp_avg(c,j,i)-dVdlog10Dp_std(c,j,i),0/))
        error_bar(c*nBins+i) = gsn_add_polyline(wks,plot,(/DpMean(i),DpMean(i)/),\
          (/top, bottom/), polyres)
        error_bar_top(c*nBins+i) = gsn_add_polyline(wks,plot,(/DpMean(i)*dfac,DpMean(i)/dfac/),\
          (/top,top/), polyres)
        error_bar_bot(c*nBins+i) = gsn_add_polyline(wks,plot,(/DpMean(i)*dfac,DpMean(i)/dfac/),\
          (/bottom,bottom/), polyres)
      end do
      
      delete(polyres@gsLineColor)
      
    end do
    draw(plot)
    frame(wks)
    
  end do
end
