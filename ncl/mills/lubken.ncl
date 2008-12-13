load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;===========================================================================

undef("createLubkenFile1999")
procedure createLubkenFile1999()

begin

  ; Read in the lubken data ...
  lubken = asciiread("analysis/Lubken1999jgr_table8.txt", (/ 59, 22 /), "float")

  ; Define the dimensions.
  z = lubken(:,0)
  z@units = "km"
  z@long_name = "Altitude"

  time = fspan(4.75, 9.75, 21)
  time@long_name = "Month"
  
  T = lubken(:, 1:21)
  T@units = "K"
  T@long_name = "Temperature"
  T!0 = "z"
  T&z = z
  T!1 = "time"
  T&time = time
  
  ; Turn it into an nc file.
  lfile = "analysis/lubken1999jgr_table8.nc"
  system("rm " + lfile)
  f = addfile(lfile, "c")
  f->time = time
  f->z = z
  f->T = T(time|:, z|:)
end

;===========================================================================

undef("createLubkenFile2004")
procedure createLubkenFile2004()

begin

  ; Read in the lubken data ...
  lubken = asciiread("analysis/Lubken2004jgr_table2.txt", (/ 59, 20 /), "float")

  ; Define the dimensions.
  z = lubken(:,0)
  z@units = "km"
  z@long_name = "Altitude"

  time = fspan(1.1, 2.9, 19)
  time@long_name = "Month"
  
  T = lubken(:, 1:19)
  T@units = "K"
  T@long_name = "Temperature"
  T!0 = "z"
  T&z = z
  T!1 = "time"
  T&time = time
  
  ; Turn it into an nc file.
  lfile = "analysis/lubken2004jgr_table2.nc"
  system("rm " + lfile)
  f = addfile(lfile, "c")
  f->time = time
  f->z = z
  f->T = T(time|:, z|:)
end


undef("plotLubkenFile1999")
procedure plotLubkenFile1999()

begin
  
  ; Open the file
  history = addfile("analysis/lubken1999jgr_table8.nc", "r")
    
  T = history->T
  
  ; Plot the data.
  wks  = gsn_open_wks("pdf", "images/lubken1999")    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  res                      = True
  res@cnFillOn             = True          ; turn on color
  res@gsnSpreadColors      = True          ; use full range of colormap
  res@lbLabelStride        = 2
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 10.0                   ; contour interval
  res@cnMinLevelValF       = 130.0             ; min level
  res@cnMaxLevelValF       = 280.0            ; max level
  res@tiMainString         = "Lubken [1999]"
    
  plot = gsn_csm_contour(wks, T(z|:, time|:), res)       ; contour the variable
end


undef("plotLubkenFile2004")
procedure plotLubkenFile2004()

begin
  
  ; Open the file
  history = addfile("analysis/lubken2004jgr_table2.nc", "r")
    
  T = history->T
  
  ; Plot the data.
  wks  = gsn_open_wks("pdf", "images/lubken2004")    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  res                      = True
  res@cnFillOn             = True          ; turn on color
  res@gsnSpreadColors      = True          ; use full range of colormap
  res@lbLabelStride        = 2
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 10.0                   ; contour interval
  res@cnMinLevelValF       = 130.0             ; min level
  res@cnMaxLevelValF       = 280.0            ; max level
  res@tiMainString         = "Lubken [2004]"
    
  plot = gsn_csm_contour(wks, T(z|:, time|:), res)       ; contour the variable
end

;===========================================================================

undef("plotTContours1999")
procedure plotTContours1999()

begin

  files = (/ "analysis/lubken1999jgr_table8.nc", \
             "analysis/125smin/tlub_1999.nc", \
             "analysis/saber_T_wg_v1.nc" \
;             "analysis/expuT125E01/tlub_1999.nc", \
;             "analysis/125expuT15/tlub_1999.nc", \
;             "analysis/125smin/tlub_1999.nc" \
;             "analysis/gw6_e125/tlub_1999.nc" \
          /)
;  names = (/ "Lubken [1999]", \
;            "tb*=1.50, eff=0.0875", \
;            "tb*=2.00, eff=0.1250", \
;            "tb*=6.00, eff=0.1250" \
;          /)
  names = (/ " Lubken [1999]", \
            " WACCM 125 levels, solar min", \
            " SABER, 2007" \
;            " WACCM 66 levels, smax", \
;            " WACCM 125 levels, smax", \
;            " WACCM 125 levels, solar min" \
          /)
          
;  timeScale = (/ 1., 30., 30., 30. /)
;  timeOffset = (/ 0., 4., 4., 1. /)
;  timeScale = (/ 1., 30., 30. /)
;  timeOffset = (/ 0., 4., 1. /)
          
  timeScale = (/ 1., 30., 30./);, 30., 30. /)
;  timeOffset = (/ 0., 1., 3., 4., 4. /)
  timeOffset = (/ 0., 4., 1. /)

  nFiles = dimsizes(files)
  plot = new(nFiles, "graphic")
  
  ; Plot the data.
  wks  = gsn_open_wks("pdf", "images/lubken1999_contours")    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  resP                      = True    
  resP@txString             = "Temperature at 70~S~o~N~N, WACCM vs. observations"
  resP@gsnMaximize          = True
  resP@gsnPanelLabelBar     = True
  resP@lbTitleString        =  "K"                ; turn on title
  resP@lbTitleOn            =  True                ; turn on title
  resP@lbTitlePosition      = "Right"              ; title position
  resP@lbTitleDirection     = "Across"             ; title direction
  resP@lbLabelStride        = 2
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@gsnPanelYWhiteSpacePercent = 5

  res                      = True
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@cnFillOn             = True          ; turn on color
  res@gsnSpreadColors      = True          ; use full range of colormap
  res@lbLabelBarOn         = False
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 10.0                   ; contour interval
  res@cnMinLevelValF       = 120.0             ; min level
  res@cnMaxLevelValF       = 270.0            ; max level
  res@cnMonoLineThickness   = False
  res@cnLineThicknesses    = (/ 1., 1., 1., 3., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. /)
;  res@trYMinF              = 35.0
;  res@trYmaxF              = 93.0
  
  ; Get the T data from each run
  print("Processing ...")
  
  do i = 0, nFiles - 1
  
    ; Open the file
    print("  " + files(i))
    history = addfile(files(i), "r")
    
    if (isStrSubset(files(i),"saber")) then
      Ttmp = dim_avg_Wrap(history->T(:, :, {70.}, :))
      z = history->lev
      timetmp = Ttmp&time
      Ttmp!1 = "z"
      Ttmp&z = z
      
      ; Temp array conforming to WACCM time dimensions, filled with missing value -999
      T = new((/180,151/), float, -999)
      T!0 = "time"
      T&time = fspan(95,274,180)
      T!1 = "z"
      T&z = z
      T({min(timetmp):max(timetmp)},:) = Ttmp(:,:)      
    else
      T = history->T
    end if
    
;    print("    dims(T):"+dimsizes(T))
    
    print("    Tmin=" + sprintf("%5.1f", min(T)))
        
    res@tiMainString       = names(i) + ", Tmin=" + sprintf("%5.2f", min(T))
    
    tmp = T(z|:, time|:)
    tmpTime = tmp&time
;    print("   tmpTime:"+min(tmpTime)+" - "+max(tmpTime))
    tmp&time = (tmpTime / timeScale(i)) + timeOffset(i)
;    print("   tmp&time:"+min(tmp&time)+" - "+max(tmp&time))
    tmp&time@long_name = "Month"
    
    plot(i) = gsn_csm_contour(wks, tmp({35.0:93.0}, {4.5:9.75}), res)       ; contour the variable

    ; Overlay the T=150 line
    ; Overlay the mesopause
    
    delete(T)
    delete(tmp)
    delete(tmpTime)
  end do

  gsn_panel(wks, plot, (/ 1, 5 /), resP)
end

;===========================================================================

undef("plotTContours2004")
procedure plotTContours2004()

begin

  files = (/ "analysis/lubken2004jgr_table2.nc", \
             "analysis/gw15_e0875/tlub_2004.nc", \
             "analysis/pmc_smin/tlub_2004.nc", \
             "analysis/gw6_e125/tlub_2004.nc" \
          /)
  names = (/ "Lubken [2004]", \
            "tb*=1.50, eff=0.0875", \
            "tb*=2.00, eff=0.1250", \
            "tb*=6.00, eff=0.1250" \
          /)
          
  timeScale = (/ 1., 30., 30., 30. /)
  timeOffset = (/ 0., -8., -8., -11.  /)
  
  nFiles = dimsizes(files)
  plot = new(nFiles, "graphic")
  
  ; Plot the data.
  wks  = gsn_open_wks("pdf", "images/lubken2004_contours")    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  resP                      = True    
  resP@txString             = "WACCM vs. Lubken [2004], 68~S~o~N~S"
  resP@gsnMaximize          = True
  resP@gsnPanelLabelBar     = True
  resP@lbTitleString        =  "K"                ; turn on title
  resP@lbTitleOn            =  True                ; turn on title
  resP@lbTitlePosition      = "Right"              ; title position
  resP@lbTitleDirection     = "Across"             ; title direction
  resP@lbLabelStride        = 2
  resP@gsnPanelXWhiteSpacePercent = 5
  resP@gsnPanelYWhiteSpacePercent = 5

  res                      = True
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@cnFillOn             = True          ; turn on color
  res@gsnSpreadColors      = True          ; use full range of colormap
  res@lbLabelBarOn         = False
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLevelSpacingF      = 10.0                   ; contour interval
  res@cnMinLevelValF       = 120.0             ; min level
  res@cnMaxLevelValF       = 270.0            ; max level
  res@cnMonoLineThickness   = False
  res@cnLineThicknesses    = (/ 1., 1., 1., 3., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. /)
;  res@trYMinF              = 35.0
;  res@trYmaxF              = 93.0
  
  ; Get the T data from each run
  print("Processing ...")
  
  do i = 0, nFiles - 1
  
    ; Open the file
    print("  " + files(i))
    history = addfile(files(i), "r")
    
    T = history->T
    
    print("    Tmin=" + sprintf("%5.1f", min(T)))
        
    res@tiMainString       = names(i) + ", Tmin=" + sprintf("%5.2f", min(T))
    
    tmp = T(z|:, time|:)
    tmpTime = tmp&time
    tmpTime@long_name = "Month"
    tmp&time = (tmpTime / timeScale(i)) + timeOffset(i)
    tmp&time@long_name = "Month"
    
    plot(i) = gsn_csm_contour(wks, tmp({35.0:93.0}, {1.1:2.9}), res)       ; contour the variable

    ; Overlay the T=150 line
    ; Overlay the mesopause
    
    delete(T)
    delete(tmp)
    delete(tmpTime)
  end do

  gsn_panel(wks, plot, (/ 1, 4 /), resP)
end

;===========================================================================

undef("plotTatAlts1999")
procedure plotTatAlts1999()

begin

  altitude = (/ 82.0, 84.0, 86.0, 88.0, 90.0, 92.0 /)

  files = (/ "analysis/lubken1999jgr_table8.nc", \
;             "analysis/expuT125E01/tlub_1999.nc", \
             "analysis/saber_T_wg_v1.nc", \
             "analysis/125expuT15/tlub_1999.nc", \
             "analysis/125smin/tlub_1999.nc" \
          /)
  names = (/ " Lubken [1999]", \
;            " WACCM 66 levels", \
            " SABER, 2007", \
            " WACCM, solar max", \
            " WACCM, solar min" \
          /)
          
;  timeScale = (/ 1., 30., 30., 30. /)
;  timeOffset = (/ 0., 4., 4., 1. /)
  timeScale = (/ 1., 30., 30., 30. /);, 30. /)
;  timeOffset = (/ 0., 3., 1., 2., 2. /)
  timeOffset = (/ 0., 1., 4., 4. /)
          
  lineColors = (/ "black", "red", "blue" , "green" /);, "orange" /)
  lineDashes = (/ 0, 0, 0 , 0 /);, 0/)
  lineThicknesses = (/ 6.0, 4.0, 4.0, 4.0 /);, 4.0 /)
  
  nFiles = dimsizes(files)
  plot = new(4, "graphic")
    
  ; Plot the data.
;  wks  = gsn_open_wks("pdf", sprintf("images/lubken1999_%3.1fkm", altitude))    ; open a ncgm file
  wks  = gsn_open_wks("pdf", "images/lubken1999_alts")    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  resP                      = True    
  resP@txString             = "WACCM vs. Observations, 70~S~o~N~N"
;  resP@gsnDraw              = False
  resP@gsnFrame             = False
  resP@gsnMaximize          = False
  resP@gsnPanelXWhiteSpacePercent = 5
;  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPaperOrientation  = "landscape"

  res                      = True
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@trXMinF              = 4.5
  res@trXMaxF              = 9.5
  res@vpYF                 = 0.85  ; Leave room for legend at the bottom.
  res@vpHeightF            = 0.7  ; Leave room for legend at the bottom.
  res@tmXBMode             = "Explicit"              ; explicit labels
  res@tmXBValues           = (/ 4., 5., 6., 7., 8., 9., 10. /) ; location of labels 
  res@tmXBLabels           = (/ 4, 5, 6, 7, 8, 9, 10 /)   ; labels themselves
  res@tiXAxisString        = "Month"
  
  resT                     = res
  resZ                     = res

  res@tiMainString       = sprintf("Temperature at %3.1f km", altitude(0))
  res@trYMinF              = 120.0
  res@trYMaxF              = 220.0

  res1                     = res
  res2                     = res
  res3                     = res
  res4                     = res
  res5                     = res

  res@tiMainString       = sprintf("Temperature at %3.1f km", altitude(0))
  res1@tiMainString       = sprintf("Temperature at %3.1f km", altitude(1))
  res2@tiMainString       = sprintf("Temperature at %3.1f km", altitude(2))
  res3@tiMainString       = sprintf("Temperature at %3.1f km", altitude(3))
  res4@tiMainString       = sprintf("Temperature at %3.1f km", altitude(4))
  res5@tiMainString       = sprintf("Temperature at %3.1f km", altitude(5))

  resT@tiMainString       = "Mesopause Temperature"
  resT@trYMinF              = 110.0
  resT@trYMaxF              = 220.0

  resZ@tiMainString       = "Mesopause Height"
  resZ@trYMinF              = 70.0
  resZ@trYMaxF              = 100.0

  ; Get the T data from each run
  print("Processing ...")
  
  do i = 0, nFiles - 1
  
    ; Open the file
    print("  " + files(i))
    history = addfile(files(i), "r")
    
    if (isStrSubset(files(i),"saber")) then
      Ttmp = dim_avg_Wrap(history->T(:, :, {70.}, :))
      z = history->lev
      timetmp = Ttmp&time
      Ttmp!1 = "z"
      Ttmp&z = z
      
      ; Temp array conforming to WACCM time dimensions, filled with missing value -999
      T = new((/180,151/), float, -999)
      T!0 = "time"
      T&time = fspan(95,274,180)
      T!1 = "z"
      T&z = z
      T({min(timetmp):max(timetmp)},:) = Ttmp(:,:)      
    else
      T = history->T
    end if
    
    print(i+" T "+dimsizes(T))
    
    print("    Tmin=" + min(T))
        
    
    tmp = T(z|:, time|:)
    tmpTime = tmp&time
    tmp@long_name = "Temperature (K)"
    tmp&time = (tmpTime / timeScale(i)) + timeOffset(i)    
    tmp&time@long_name = "Month"
    
    ; Determine the location of the temperature minimum at each
    ; time.
    nTimes = dimsizes(tmp&time)
    nZs = dimsizes(tmp&z)
    imp = new(nTimes, "integer")
    dataT = new(nTimes, "float")
    Zmin = new(nTimes, "float")
    
    print(i+" nZs "+nZs)
        
    do j = 0, nTimes-1
      imp(j) = nZs - 1 - minind(tmp(::-1, j))
      if ismissing(imp(j)) then
        dataT(j) = -999
        Zmin(j) = -999
      else
        dataT(j) = tmp(imp(j), j)
        Zmin(j) = tmp&z(imp(j))
      end if
    end do
    
    res@xyLineColor = lineColors(i)
    res1@xyLineColor = lineColors(i)
    res2@xyLineColor = lineColors(i)
    res3@xyLineColor = lineColors(i)
    res4@xyLineColor = lineColors(i)
    res5@xyLineColor = lineColors(i)
    resT@xyLineColor = lineColors(i)
    resZ@xyLineColor = lineColors(i)
    
    res@xyDashPatterns = lineDashes(i)
    res1@xyDashPatterns = lineDashes(i)
    res2@xyDashPatterns = lineDashes(i)
    res3@xyDashPatterns = lineDashes(i)
    res4@xyDashPatterns = lineDashes(i)
    res5@xyDashPatterns = lineDashes(i)
    resT@xyDashPatterns = lineDashes(i)
    resZ@xyDashPatterns = lineDashes(i)
    
    res@xyLineThicknessF = lineThicknesses(i)
    res1@xyLineThicknessF = lineThicknesses(i)
    res2@xyLineThicknessF = lineThicknesses(i)
    res3@xyLineThicknessF = lineThicknesses(i)
    res4@xyLineThicknessF = lineThicknesses(i)
    res5@xyLineThicknessF = lineThicknesses(i)
    resT@xyLineThicknessF = lineThicknesses(i)
    resZ@xyLineThicknessF = lineThicknesses(i)

    if (i .eq. 0) then
      plot(0) = gsn_csm_xy(wks, tmp&time(:), tmp({altitude(0)}, :), res)
      plot(1) = gsn_csm_xy(wks, tmp&time(:), tmp({altitude(1)}, :), res1)
      plot(2) = gsn_csm_xy(wks, tmp&time(:), tmp({altitude(2)}, :), res2)
;      plot(3) = gsn_csm_xy(wks, tmp&time(:), tmp({altitude(3)}, :), res3)
;      plot(4) = gsn_csm_xy(wks, tmp&time(:), tmp({altitude(4)}, :), res4)
;      plot(5) = gsn_csm_xy(wks, tmp&time(:), tmp({altitude(5)}, :), res5)
      plot(3) = gsn_csm_xy(wks, tmp&time(:), dataT(:), resT)
;      plot(7) = gsn_csm_xy(wks, tmp&time(:), tmp&z(imp), resZ)
    else
      overlay(plot(0), gsn_csm_xy(wks, tmp&time(:), tmp({altitude(0)}, :), res))
      overlay(plot(1), gsn_csm_xy(wks, tmp&time(:), tmp({altitude(1)}, :), res))
      overlay(plot(2), gsn_csm_xy(wks, tmp&time(:), tmp({altitude(2)}, :), res))
;      overlay(plot(3), gsn_csm_xy(wks, tmp&time(:), tmp({altitude(3)}, :), res))
;      overlay(plot(4), gsn_csm_xy(wks, tmp&time(:), tmp({altitude(4)}, :), res))
;      overlay(plot(5), gsn_csm_xy(wks, tmp&time(:), tmp({altitude(5)}, :), res))
      overlay(plot(3), gsn_csm_xy(wks, tmp&time(:), dataT(:), resT))
;      overlay(plot(7), gsn_csm_xy(wks, tmp&time(:), Zmin(:), resZ))
    end if

    delete(T)
    delete(tmp)
    delete(tmpTime)
    delete(imp)
    delete(dataT)
    delete(Zmin)
  end do
  
  gsn_panel(wks, plot, (/ 2, 2 /), resP)

  legend = create "Legend" legendClass wks 
    "vpXF"                     : 0.1                    ; orientation on page
    "vpYF"                     : -0.275
    "vpWidthF"                 : 0.8                     ; width
    "vpHeightF"                : 0.1                     ; height
    "lgPerimOn"                : True                   ; no perimeter
    "lgItemCount"              : nFiles                      ; how many
    "lgLineLabelStrings"       : names
    "lgLabelsOn"               : False                   ; no default lables
    "lgMonoLineThickness"      : False
    "lgLineThicknesses"        : (/ 6.0, 4.0, 4.0, 4.0, 4.0, 4.0 /)             ; line paters
    "lgDashIndexes"            : lineDashes             ; line paters
    "lgLineColors"             : lineColors
    "lgMonoLineLabelFontColor" : True                    ; one label color
    "lgAutoManage"             : False
    "lgLabelPosition"          : "Center"
    "lgLabelJust"              : "CenterCenter"
  end create

  draw(legend)  
  frame(wks)

end

;===========================================================================

undef("plotTatZ1999")
procedure plotTatZ1999(altitude:float)

begin

  files = (/ "analysis/lubken1999jgr_table8.nc", \
             "analysis/saber_T_wg_v1.nc", \
             "analysis/expuT125E01/tlub_1999.nc" \
;             "analysis/pmc_smin/tlub_1999.nc" \
;             "analysis/gw6_e125/tlub_1999.nc" \
          /)
  names = (/ " Lubken [1999]", \
            " SABER, 2007", \
            " WACCM 66 levels" \
;            "Default (tb*=6.00, eff=0.1250)" \
          /)
          
;  timeScale = (/ 1., 30., 30., 30. /)
;  timeOffset = (/ 0., 4., 4., 1. /)
  timeScale = (/ 1., 30., 30. /)
  timeOffset = (/ 0., 3., 3. /)
          
  lineColors = (/ "black", "red", "blue"/) ;, "green" /)
  lineDashes = (/ 0, 0, 0 /)
  lineThicknesses = (/ 6.0, 4.0, 4.0 /)
  
  nFiles = dimsizes(files)
  plot = new(3, "graphic")
    
  ; Plot the data.
  wks  = gsn_open_wks("pdf", sprintf("images/lubken1999_%3.1fkm", altitude))    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  resP                      = True    
  resP@txString             = "WACCM vs. Lubken [1999], 70~S~o~N~N"
;  resP@gsnDraw              = False
  resP@gsnFrame             = False
  resP@gsnMaximize          = True
  resP@gsnPanelXWhiteSpacePercent = 5
;  resP@gsnPanelYWhiteSpacePercent = 5

  res                      = True
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@trXMinF              = 4.5
  res@trXMaxF              = 9.5
  res@vpYF                 = 0.85  ; Leave room for legend at the bottom.
  res@vpHeightF            = 0.7  ; Leave room for legend at the bottom.
  res@tmXBMode             = "Explicit"              ; explicit labels
  res@tmXBValues           = (/ 4., 5., 6., 7., 8., 9., 10. /) ; location of labels 
  res@tmXBLabels           = (/ 4, 5, 6, 7, 8, 9, 10 /)   ; labels themselves
  res@tiXAxisString        = "Month"
  
  resT                     = res
  resZ                     = res


  res@tiMainString       = sprintf("Temperature at %3.1f km", altitude)
  res@trYMinF              = 120.0
  res@trYMaxF              = 220.0


  resT@tiMainString       = "Mesopause Temperature"
  resT@trYMinF              = 110.0
  resT@trYMaxF              = 220.0

  resZ@tiMainString       = "Mesopause Height"
  resZ@trYMinF              = 70.0
  resZ@trYMaxF              = 100.0

  ; Get the T data from each run
  print("Processing ...")
  
  do i = 0, nFiles - 1
  
    ; Open the file
    print("  " + files(i))
    history = addfile(files(i), "r")
    
    T = history->T
    
    print("    Tmin=" + min(T))
        
    
    tmp = T(z|:, time|:)
    tmpTime = tmp&time
    tmp@long_name = "Temperature (K)"
    tmp&time = (tmpTime / timeScale(i)) + timeOffset(i)    
    tmp&time@long_name = "Month"
    
    ; Determine the location of the temperature minimum at each
    ; time.
    nTimes = dimsizes(tmp&time)
    nZs = dimsizes(tmp&z)
    imp = new(nTimes, "integer")
    dataT = new(nTimes, "float")
    
    do j = 0, nTimes-1
      imp(j) = nZs - 1 - minind(tmp(::-1, j))
      dataT(j) = tmp(imp(j), j)
    end do
    
    res@xyLineColor = lineColors(i)
    resT@xyLineColor = lineColors(i)
    resZ@xyLineColor = lineColors(i)
    
    res@xyDashPatterns = lineDashes(i)
    resT@xyDashPatterns = lineDashes(i)
    resZ@xyDashPatterns = lineDashes(i)
    
    res@xyLineThicknessF = lineThicknesses(i)
    resT@xyLineThicknessF = lineThicknesses(i)
    resZ@xyLineThicknessF = lineThicknesses(i)

    if (i .eq. 0) then
      plot(0) = gsn_csm_xy(wks, tmp&time(:), tmp({altitude}, :), res)
      plot(1) = gsn_csm_xy(wks, tmp&time(:), dataT(:), resT)
      plot(2) = gsn_csm_xy(wks, tmp&time(:), tmp&z(imp), resZ)
    else
      overlay(plot(0), gsn_csm_xy(wks, tmp&time(:), tmp({altitude}, :), res))
      overlay(plot(1), gsn_csm_xy(wks, tmp&time(:), dataT(:), resT))
      overlay(plot(2), gsn_csm_xy(wks, tmp&time(:), tmp&z(imp), resZ))
    end if

    delete(T)
    delete(tmp)
    delete(tmpTime)
    delete(imp)
    delete(dataT)
  end do
  
  gsn_panel(wks, plot, (/ 1, 3 /), resP)

  legend = create "Legend" legendClass wks 
    "vpXF"                     : 0.1                    ; orientation on page
    "vpYF"                     : 0.275
    "vpWidthF"                 : 0.8                     ; width
    "vpHeightF"                : 0.1                     ; height
    "lgPerimOn"                : True                   ; no perimeter
    "lgItemCount"              : nFiles                      ; how many
    "lgLineLabelStrings"       : names
    "lgLabelsOn"               : False                   ; no default lables
    "lgMonoLineThickness"      : False
    "lgLineThicknesses"        : (/ 6.0, 4.0, 4.0, 4.0, 4.0, 4.0 /)             ; line paters
    "lgDashIndexes"            : lineDashes             ; line paters
    "lgLineColors"             : lineColors
    "lgMonoLineLabelFontColor" : True                    ; one label color
    "lgAutoManage"             : False
    "lgLabelPosition"          : "Center"
    "lgLabelJust"              : "CenterCenter"
  end create

  draw(legend)  
  frame(wks)

end

;===========================================================================

undef("plotTatZ2004")
procedure plotTatZ2004(altitude:float)

begin

  files = (/ "analysis/lubken2004jgr_table2.nc", \
             "analysis/gw15_e0875/tlub_2004.nc", \
             "analysis/pmc_smin/tlub_2004.nc", \
             "analysis/gw6_e125/tlub_2004.nc" \
          /)
  names = (/ "Lubken [2004]", \
            "tb*=1.50, eff=0.0875", \
            "tb*=2.00, eff=0.1250", \
            "tb*=6.00, eff=0.1250" \
          /)
          
  timeScale = (/ 1., 30., 30., 30. /)
  timeOffset = (/ 0., -8., -8., -11.  /)
          
  lineColors = (/ "black", "blue", "green", "red" /)
  lineDashes = (/ 0, 0, 0, 0 /)
  lineThicknesses = (/ 6.0, 4.0, 4.0, 4.0 /)
  
  nFiles = dimsizes(files)
  plot = new(3, "graphic")
    
  ; Plot the data.
  wks  = gsn_open_wks("pdf", "images/lubken2004_82km")    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow")  ; choose a colormap

  resP                      = True    
  resP@txString             = "WACCM vs. Lubken [2004], 68~S~o~N~S"
  resP@gsnFrame             = False
  resP@gsnMaximize          = True
  resP@gsnPanelXWhiteSpacePercent = 5
;  resP@gsnPanelYWhiteSpacePercent = 5

  res                      = True
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@trXMinF              = -1.0
  res@trXMaxF              = 3.0
  res@vpYF                 = 0.85  ; Leave room for legend at the bottom.
  res@vpHeightF            = 0.7  ; Leave room for legend at the bottom.
  res@tmXBMode             = "Explicit"              ; explicit labels
  res@tmXBValues           = (/ -1., 0., 1., 2., 3. /) ; location of labels 
  res@tmXBLabels           = (/ 11, 12, 1, 2, 3 /)   ; labels themselves
  res@tiXAxisString        = "Month"
  
  resT                     = res
  resZ                     = res


  res@tiMainString       = "Temperature at 82 km"
  res@trYMinF              = 120.0
  res@trYMaxF              = 220.0


  resT@tiMainString       = "Mesopause Temperature"
  resT@trYMinF              = 110.0
  resT@trYMaxF              = 220.0

  resZ@tiMainString       = "Mesopause Height"
  resZ@trYMinF              = 70.0
  resZ@trYMaxF              = 100.0

  ; Get the T data from each run
  print("Processing ...")
  
  do i = 0, nFiles - 1
  
    ; Open the file
    print("  " + files(i))
    history = addfile(files(i), "r")
    
    T = history->T
    
    print("    time(0)=" + T&time(0))
    print("    Tmin=" + min(T))
        
    
    tmp = T(z|:, time|:)
    tmpTime = tmp&time
    tmp@long_name = "Temperature (K)"
    tmpTime@long_name = "Month"
    tmp&time = (tmpTime / timeScale(i)) + timeOffset(i)
    
    print(tmp&time)
    
    tmp&time@long_name = "Month"
    
    ; Determine the location of the temperature minimum at each
    ; time.
    nTimes = dimsizes(tmp&time)
    nZs = dimsizes(tmp&z)
    imp = new(nTimes, "integer")
    dataT = new(nTimes, "float")
    
    do j = 0, nTimes-1
      imp(j) = nZs - 1 - minind(tmp(::-1, j))
      dataT(j) = tmp(imp(j), j)
    end do
    
;    print(imp)
;    print(tmp&z(imp))

    res@xyLineColor = lineColors(i)
    resT@xyLineColor = lineColors(i)
    resZ@xyLineColor = lineColors(i)
    
    res@xyDashPatterns = lineDashes(i)
    resT@xyDashPatterns = lineDashes(i)
    resZ@xyDashPatterns = lineDashes(i)
    
    res@xyLineThicknessF = lineThicknesses(i)
    resT@xyLineThicknessF = lineThicknesses(i)
    resZ@xyLineThicknessF = lineThicknesses(i)

    if (i .eq. 0) then
      plot(0) = gsn_csm_xy(wks, tmp&time(:), tmp({altitude}, :), res)
      plot(1) = gsn_csm_xy(wks, tmp&time(:), dataT(:), resT)
      plot(2) = gsn_csm_xy(wks, tmp&time(:), tmp&z(imp), resZ)
    else
      overlay(plot(0), gsn_csm_xy(wks, tmp&time(:), tmp({altitude}, :), res))
      overlay(plot(1), gsn_csm_xy(wks, tmp&time(:), dataT(:), resT))
      overlay(plot(2), gsn_csm_xy(wks, tmp&time(:), tmp&z(imp), resZ))
    end if

    delete(T)
    delete(tmp)
    delete(tmpTime)
    delete(imp)
    delete(dataT)
  end do
  
  gsn_panel(wks, plot, (/ 1, 3 /), resP)

  legend = create "Legend" legendClass wks 
    "vpXF"                     : 0.1                    ; orientation on page
    "vpYF"                     : 0.275
    "vpWidthF"                 : 0.8                    ; width
    "vpHeightF"                : 0.1                  ; height
    "lgPerimOn"                : True                   ; no perimeter
    "lgItemCount"              : nFiles                 ; how many
    "lgLineLabelStrings"       : names
    "lgLabelsOn"               : False                   ; no default lables
    "lgMonoLineThickness"      : False
    "lgLineThicknesses"        : lineThicknesses             ; line paters
    "lgDashIndexes"            : lineDashes             ; line paters
    "lgLineColors"             : lineColors
    "lgMonoLineLabelFontColor" : True                    ; one label color
    "lgAutoManage"             : False
    "lgLabelPosition"          : "Center"
    "lgLabelJust"              : "CenterCenter"
  end create

  draw(legend)  
  frame(wks)

end

;===========================================================================

;createLubkenFile1999()
;createLubkenFile2004()

;plotLubkenFile1999()
;plotLubkenFile2004()
;plotTContours1999()
;plotTContours2004()
;plotTatZ1999(82.0)
;plotTatZ2004(82.0)
;plotTatZ1999(86.0)
;plotTatZ1999(88.0)
;plotTatZ1999(90.0)
;plotTatZ1999(92.0)
plotTatAlts1999()

