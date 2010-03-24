load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotvgSootCol")
procedure plotvgSootCol(SOOT01_col[*][*]:float) 
local O3_col_diff,wks,res,lev
begin 

  lat=SOOT01_col&lat
  
  xdim=dimsizes(SOOT01_col)
  print("xdim="+xdim)
  ntimes=min((/xdim(0)/))
    
  wks = gsn_open_wks("pdf", "vgSootColumnLatSeason")
;  gsn_define_colormap(wks,"BlRe")
;  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  gsn_define_colormap(wks,"rainbow")
  res                     = True          ; plot mods desired

  res@cnFillOn            = True          ; turn on color
;  res@lbLabelBarOn        = False         ; turn off color labelbar
  res@gsnSpreadColors     = True          ; use full range of colormap
  res@gsnMaximize         = True          ; expand plot size to maximum
  res@gsnPaperOrientation = "landscape"
  res@vpWidthF            = 1.2           ; change aspect ratio of plot
  res@vpHeightF           = 0.7
;  res@gsnCenterStringFontHeightF = 24
;  res@gsnRightStringFontHeightF = 24
  res@tiYAxisString       = "Latitude"
  res@cnLineLabelPlacementMode= "Computed"
  res@cnLineLabelDensityF = 5
  res@tmXBLabelFontHeightF =0.03
  res@tmYLLabelFontHeightF =0.03
  
  res@tmXBLabels =(/"J","F","M","A","M","J","J","A","S","O","N","D"/)   
  res@tmXBMode           = "Explicit"         ; Define own tick mark labels.
  res@tmXBValues         = SOOT01_col&time(1:12)
  
  res@tmXBMinorValues    = SOOT01_col&time(0:12)+0.5
  res@tmXBMinorOutwardLengthF = 0.02    ; Increase length and size of
  res@tmXBMinorLengthF        = 0.02    ; minor tickmarks
  res@tmXBMinorThicknessF     = 2.0
;  res@tmXBMajorThicknessF     = 2.0
;  res@tmXBMajorOutwardLengthF = 0.01
  res@tmXBMajorLengthF        = 0.00
;===================================================================  
  res@trXMinF = 0.5   
  res@trXMaxF = 12.5  
; Define our own contour levels.  
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnMinLevelValF=400
  res@cnMaxLevelValF=4400
  res@cnLevelSpacingF     = 400.
  res@gsnCenterString      = "Soot Column (~F8~m~F21~g m~S~-2~N~)"
  res@gsnRightString      = " "
  res@lbLabelAngleF       = 60
  plot = gsn_csm_lat_time(wks,SOOT01_col(lat|:, time|:),res)      ; contour the variable
  
end


undef("readvgSootCol")
procedure readvgSootCol(rootPath:string,historyRun:string)

begin
  print((/historyRun/))
  historyPattern = rootPath+"/"+historyRun+"/h0zm/"+historyRun+".cam2.h0zm.[01][0-9].nc"
  historyNames = systemfunc("ls " + historyPattern)
;  print((/historyNames/))
    
  a=addfile(historyNames(0),"r")
  SOOT01_col=new((/14,dimsizes(a->lat)/),float)
  SOOT01_col!0 = "time"
  SOOT01_col!1 = "lat"
  SOOT01_col&lat = a->lat
  print(dimsizes(SOOT01_col))
  
  time = new(14,integer)
  
  do i = 1,12
    time(i)=i
    a=addfile(historyNames(i-1),"r")
    SOOT01_col(i,:)=a->SOOT01_col
  end do
  time(13)=13
  time(0)=0
  
  SOOT01_col(0,:)=SOOT01_col(12,:)
  SOOT01_col(13,:)=SOOT01_col(1,:)
  
  SOOT01_col&time=time
  
  if (SOOT01_col@units .eq. "g/cm2") then
    SOOT01_col@units="ug m~-2"
    SOOT01_col=SOOT01_col*1e10
  end if

  print("plotting ozone column...")
  plotvgSootCol(SOOT01_col)
  
end

historyRun="vgsoot08"
;rootPath="/data5/mmills"
rootPath="/Volumes/Data/Models/waccm/carma-soot/run"

readvgSootCol(rootPath,historyRun)
