load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotvgSootCol")
procedure plotvgSootCol(SOOTMD_col[*][*]:float) 
local O3_col_diff,wks,res,lev
begin 

  lat=SOOTMD_col&lat
  
  xdim=dimsizes(SOOTMD_col)
  print("xdim="+xdim)
  ntimes=min((/xdim(0)/))
    
  wks = gsn_open_wks("pdf", "$IMAGE_PATH/vgSootColumnLatSeason")
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
  res@tmXBValues         = SOOTMD_col&time(1:12)
  
  res@tmXBMinorValues    = SOOTMD_col&time(0:12)+0.5
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
  res@cnMinLevelValF=750
  res@cnMaxLevelValF=9000
  res@cnLevelSpacingF     = 750.
  res@gsnCenterString      = "Soot Column (ng m~S~-2~N~)"
  res@gsnRightString      = " "
  res@lbLabelAngleF       = 60
  plot = gsn_csm_lat_time(wks,SOOTMD_col(lat|:, time|:),res)      ; contour the variable
  
end


undef("readvgSootCol")
procedure readvgSootCol(rootPath:string,historyRun:string)

begin
  print((/historyRun/))
  historyPattern = rootPath+"/"+historyRun+"/h0zm/"+historyRun+".cam2.h0zm.[01][0-9].nc"
  historyNames = systemfunc("ls " + historyPattern)
;  print((/historyNames/))
    
  a=addfile(historyNames(0),"r")
  gw=a->gw ; Gaussian weights (latitude)
  
  SOOTMD_col=new((/14,dimsizes(a->lat)/),float)
  SOOTMD_col!0 = "time"
  SOOTMD_col!1 = "lat"
  SOOTMD_col&lat = a->lat
  print(dimsizes(SOOTMD_col))
  
  SOOT_mass=new(12,double)
  SOOT_mass_SH=new(12,double)
    
  time = new(14,integer)
  
  do i = 1,12
    time(i)=i
    a=addfile(historyNames(i-1),"r")
    sootcol=a->SOOTMD_col
    SOOTMD_col(i,:)=sootcol
    sootmass=sootcol(:)#gw(:)
    sootmassSH=sootcol({-90:-1})#gw({-90:-1})
    fracSH=sootmassSH/sootmass
    
    sootmass=sootmass*510065600*1e6 ; Surface area of the Earth (m2)
    sootmass=sootmass/1e9 ; ng->g
    sootmass=sootmass/1e6 ; g->tonnes
    
    sootmassSH=sootmassSH*510065600*1e6 ; Surface area of the Earth (m2)
    sootmassSH=sootmassSH/1e9 ; ng->g
    sootmassSH=sootmassSH/1e6 ; g->tonnes
    
    print(i+" "+(/sootmass/)+" "+(/sootmassSH/)+" "+(/fracSH/))
;    print((/sootmass/))
    SOOT_mass(i-1)=(/sootmass/)
    SOOT_mass_SH(i-1)=(/sootmassSH/)
  end do
  
  print("mean soot mass:"+avg(SOOT_mass))
  print("SH soot mass:"+avg(SOOT_mass_SH))
  delete(fracSH)
  fracSH=avg(SOOT_mass_SH)/avg(SOOT_mass)
  print("SH frac:"+(/fracSH/))
  
  time(13)=13
  time(0)=0
  
  SOOTMD_col(0,:)=SOOTMD_col(12,:)
  SOOTMD_col(13,:)=SOOTMD_col(1,:)
  
  SOOTMD_col&time=time
  
;  if (SOOTMD_col@units .eq. "g/cm2") then
;    SOOTMD_col@units="ng m~-2"
;    SOOTMD_col=SOOTMD_col*1e13
;  end if

;  print("plotting soot column...")
;  plotvgSootCol(SOOTMD_col)
  
end

historyRun="vgsoot08"
;rootPath="/data5/mmills"
rootPath="/Volumes/Data/Models/waccm/run/carma-soot"

readvgSootCol(rootPath,historyRun)
