load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************

begin

cases = (/"b.e10.BRCP45WCN.f19_g16.nsplit008.001", \
          "b.e10.BRCP45WCN.f19_g16.nsplit016.001", \          
          "b.e10.BRCP45WCN.f19_g16.nsplit032.001", \
          "b.e10.BRCP45WCN.f19_g16.nsplit064.001", \          
          "b.e10.BRCP45WCN.f19_g16.nsplit128.001" /)
          
rootpath="/glade/scratch/mmills/archive/"

ncases = dimsizes(cases)

do i = 0, ncases-1
  filename = rootpath + cases(i) +"/atm/proc/h0gm/"+cases(i)+".cam.h0.globavg.nc"
  print("reading "+filename)
  fileh = addfile(filename,"r")
  if (i.eq.0) then
    time=fileh->time
    ntimes=dimsizes(time)
    FSNS=new((/ncases,ntimes/),double)
    PRECT=FSNS
    CLDLOW=FSNS
    CLDMED=FSNS
    CLDHGH=FSNS
    CLDTOT=FSNS
  end if
  FSNS(i,:)=fileh->FSNS
  PRECT(i,:)=fileh->PRECT
  CLDLOW(i,:)=fileh->CLDLOW
  CLDMED(i,:)=fileh->CLDMED
  CLDHGH(i,:)=fileh->CLDHGH
  CLDTOT(i,:)=fileh->CLDTOT
  delete(fileh)
  delete(filename)
end do      

;************************************************
; create plots
;************************************************
  
  wks = gsn_open_wks("pdf","$IMG_PATH/nsplit")               ; open a ps file
  plot = new(6,graphic)                          ; create a plot array
  
  res          = True
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@xyMonoLineColor   = False
  res@xyLineColors      = (/"black","purple","blue","green","red"/)
  res@xyMonoDashPattern   = True
  res@xyLineThicknessF = 4.0
;  res@trXmaxF=31.0

  plot(0)=gsn_csm_xy(wks,time,FSNS,res) 
  plot(1)=gsn_csm_xy(wks,time,PRECT,res)
  plot(2)=gsn_csm_xy(wks,time,CLDLOW,res) 
  plot(3)=gsn_csm_xy(wks,time,CLDMED,res)
  plot(4)=gsn_csm_xy(wks,time,CLDHGH,res) 
  plot(5)=gsn_csm_xy(wks,time,CLDTOT,res) 
;***********************************************
; legend resources
;***********************************************
  lgres                    = True
  lgres@lgLineColors     = res@xyLineColors
  lgres@lgLineThicknessF   = res@xyLineThicknessF
  lgres@lgMonoDashIndex = True
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; set the legend label font thickness
  lgres@vpWidthF           = 0.2           ; width of legend (NDC)
  lgres@vpHeightF          = 0.2            ; height of legend (NDC)
;  lgres@lgPerimColor       = "orange"       ; draw the box perimeter in orange
;  lgres@lgPerimThicknessF  = 5.0            ; thicken the box perimeter
  lbid = gsn_create_legend(wks,5,(/"8","16","32","64","128"/),lgres)         ; create legend

  amres = True
  amres@amParallelPosF   = 0.25	; move legend to the left
  amres@amOrthogonalPosF = 0.3          ; move the legend down
  annoid1 = gsn_add_annotation(plot(1),lbid,amres)   ; add legend to plot
  
;************************************************
; create panel
;************************************************
; draw panel with white space added
 resP                 = True
 resP@gsnPanelYWhiteSpacePercent = 5
 resP@gsnPanelXWhiteSpacePercent = 10
 resP@gsnMaximize = True
 gsn_panel(wks,plot,(/3,2/),resP)

end  
  
  
