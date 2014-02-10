load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  rootpath = "/glade/scratch/mmills/archive/"
  
  caseid = (/  \
              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.QBOon.003", \
              "b.e12.B55TRW5CN.f19_g16.Pin10Tg.20-25km.15S-15N.zavg.QBOon.003", \
              "b.e12.BRCP45W5CN.f19_g16.MAMstratchem.003" /)
  
  yearString = (/ "1992-93", "1996-97", "2006-07" /)
  
  nCases = dimsizes(caseid)
  
  do i = 0, nCases-1
  
    a=addfile(rootpath+caseid(i)+"/atm/proc/SizeDistN2O_"+yearString(i)+".nc","r")
    if (i.eq.0) then
      DpMean=a->diameter
      nBins=dimsizes(DpMean)
      N2Obounds=a->N2Oppb
      nRanges=dimsizes(N2Obounds)-1
      dVdlog10Dp_avg=new((/nCases,nRanges,nBins/),"float")
      dVdlog10Dp_avg!2="diameter"
      dVdlog10Dp_avg&diameter=DpMean
      dVdlog10Dp_std=dVdlog10Dp_avg
    end if
    
    dVdlog10Dp_avg(i,:,:)=a->dVdlog10Dp_avg
    dVdlog10Dp_std(i,:,:)=a->dVdlog10Dp_std
  
  end do
              
  pdffile="SizeDistN2O_multi"
  wks   = gsn_open_wks ("pdf", "$IMAGE_PATH/"+pdffile)                   ; open workstation
  res=True
  res@gsnFrame               = False                     ; don't draw yet
  res@gsnDraw                = False                     ; don't advance frame
  res@xyXStyle="Log"
  res@trXMinF=0.05
  res@trXMaxF=5.0
  res@trYMinF=-0.1
  res@trYMaxF=3.0
  res@xyMonoLineColor = False
  res@xyLineColors    = (/"red","green","blue"/)
  res@xyMonoDashPattern = True
  
  polyres                   = True                       ; marker resources
  polyres@gsMarkerIndex     = 1                          ; index for circle
  polyres@gsMarkerSizeF     = .02                        ; size 
  
  do j = 0, nRanges-1
  
    plot = gsn_csm_xy (wks, DpMean, dVdlog10Dp_avg(:,j,:), res)

  ;************************************************
  ; add error bars
  ;************************************************ 
    error_bar = new(nCases*nBins,graphic)
    do c = 0, nCases-1 
      polyres@gsLineColor = res@xyLineColors(c)

      do i = 0, nBins-1
        error_bar(c*nBins+i) = gsn_add_polyline(wks,plot,(/DpMean(i),DpMean(i)/),\
          (/dVdlog10Dp_avg(c,j,i)+dVdlog10Dp_std(c,j,i),\
          max((/dVdlog10Dp_avg(c,j,i)-dVdlog10Dp_std(c,j,i),0/))/),\
          polyres)
      end do
      
      delete(polyres@gsLineColor)
      
    end do
    draw(plot)
    frame(wks)
    
  end do
end
