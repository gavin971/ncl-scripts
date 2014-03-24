load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  rootpath = "/glade/scratch/mmills/archive/"
  
;  caseid = (/  \
;              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.QBOon.003", \
;;              "b.e12.B55TRW5CN.f19_g16.Pin15Tg.16-20km.QBOon.001", \
;              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.QBOon.003", \
;              "b.e12.BRCP45W5CN.f19_g16.MAMstratchem.003" /)
;  caseid = (/  \
;              "f.e12.FSDW5.f19_f19.Pin15Tg.20-27km.0-18N.96-118E.QBOon.001", \
;              "f.e12.FSDW5.f19_f19.Pin15Tg.20-27km.0-18N.96-118E.QBOon.001", \
;              "b.e12.BRCP45W5CN.f19_g16.MAMstratchem.003" /)
;  caseid = (/  \
;              "f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.QBOon.002", \
;              "f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.QBOon.002", \
;              "b.e12.BRCP45W5CN.f19_g16.MAMstratchem.003" /)
;  caseid = (/  \
;              "f.e12.FSDW5.f19_f19.Pin20Tg.19-27km.0-20N.120E.QBOon.001", \
;              "f.e12.FSDW5.f19_f19.Pin20Tg.19-27km.0-20N.120E.QBOon.001", \
;              "b.e12.BRCP45W5CN.f19_g16.MAMstratchem.003" /)
  caseid = (/  \
              "f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.a3dgnl44.001", \
              "f.e12.FSDW5.f19_f19.Pin15Tg.19-27km.0N-20N.120E.a3dgnl44.001", \
              "f.e12.FSDW5.f19_f19.Pin10Tg.20-25km.15S-15N.zavg.a3dgnl44.001", \
              "f.e12.FSDW5.f19_f19.Pin15Tg.19-27km.0N-20N.120E.a3dgnl44.001" \
               /)
  
;  yearString = (/ "1991_06-12", "1992", "2006-07" /)
  yearString = (/ "1992-93", "1992-93", "1996-97", "1996-97" /)
  
  nCases = dimsizes(caseid)
  
  do i = 0, nCases-1
  
    a=addfile(rootpath+caseid(i)+"/atm/proc/SizeDistN2O_"+yearString(i)+".nc","r")
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
              
  pdffile="SizeDistN2O_multi_"+caseid(0)+"_c"+datestr
  wks   = gsn_open_wks ("pdf", "$IMAGE_PATH/"+pdffile)                   ; open workstation
  res=True
  res@vpHeightF             = 0.4              ; change aspect ratio
  res@vpWidthF              = 0.7
  res@gsnFrame               = False                     ; don't draw yet
  res@gsnDraw                = False                     ; don't advance frame
  res@xyXStyle="Log"
  res@trXMinF=0.05
  res@trXMaxF=5.0
  res@trYMinF=-0.2
  res@trYMaxF=3.1
  res@xyMonoLineColor = False
  res@xyLineColors    = (/"red","red","green3","green3","blue"/)
  res@xyMonoDashPattern = False
  res@xyDashPatterns = (/1,0,1,0,0/)
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
