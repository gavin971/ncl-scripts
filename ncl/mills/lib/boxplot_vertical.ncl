; ----------------------------------------------------------
; Mike Mills, adapted from boxplot by Adam Phillips
;
undef("boxplot_vertical")
function boxplot_vertical (wks:graphic,x[*][*][*]:numeric,y[*]:numeric,boxOpts:logical,plotres:logical,lineres:logical)
;
;    This function creates a blank linLog plot object, on which box plots are created by extensive
;    use of gsn_add_polyline. The user can draw as many or as few box plots as they wish.
;    The plot is not drawn, and the frame is not advanced. May be used with gsn_panel. Missing data
; may be present in the input Y array, but not the input X array.
;    The function options are as follows:
;
;    wks     ----- The workstation
;    x[*]    ----- A one-dimensional array containing the X-axis values of where the box plot(s)
;                  will be drawn.
;    y[*][*] ----- A two-dimensional array, where the rightmost dimension contains the box plot
;               reference pt. data. y(n,0)=bottom_value, y(n,1)=bottom_value_of_box,
;               y(n,2)=mid-value_of_box,y(n,3)=top_value_of_box,y(n,4)=top_value
;    boxOpts ----- 2 options attached as attributes may be attached here.
;               boxWidth ---- Scalar or array that contains the widths of the boxes.
;               boxColors ---- Scalar or array that contains the colors that the boxes will be drawn in.
;                           Ex. boxOpts@boxColors =  (/"green","blue"/)
;                        If the number of colors specified does not equal the number of
;                        boxes drawn, only the first color specified will be used.
;    plotres ----- An optional xy-plot resource list. Will override defaults where applicable.
;    lineres ----- An optional resource list pertaining to the lines used to draw the boxes/lines.
;
;    Example call: plot3 = boxplot_vertical(wks,ntime,newyval,opti,res,False)
;               draw(plot3)
;               frame(wks)
;
begin
     dimquery = dimsizes(x)
     numseries = dimquery(0)
     numbox = dimquery(1)
     print("boxplot_vertical: numbox="+numbox+" numseries="+numseries)
     boxWidths = new((/numbox/),float)
     if (numbox.ne.dimsizes(y)) then
        print("boxplot_vertical: Fatal: Y must be one-dimensional must have the same size as X dimension 1")
        exit
     end if
     if (any(ismissing(y))) then
        print("boxplot_vertical: Fatal: Y array cannot contain missing data, exiting")
        exit
     end if

;    Developing y-axis
     yAxis = new(numbox+2,typeof(y))
     yAxis(1:numbox) = y
     if (numbox.ne.1) then
        dy = y(1)-y(0)
        yAxis(0) = y(0)-dy
        yAxis(numbox+1) = y(numbox-1)+dy
     else
        dy = 1
        yAxis(0) = y-dy
        yAxis(2) = y+dy
     end if
;     print(yAxis)

     if (boxOpts) then
        if (isatt(boxOpts,"boxWidth")) then
           if (dimsizes(boxOpts@boxWidth).ne.1.and.dimsizes(boxOpts@boxWidth).ne.numbox) then
              print("boxplot_vertical: Number of input box widths must either equal 1 or the number of boxes ("+numbox+"). Using first specified box width only.")
              boxWidths(:) = boxOpts@boxWidth(0)
           else
              boxWidths = boxOpts@boxWidth
           end if
        else
           boxWidths(:) = dy*.3
        end if
     else
        boxWidths(:) = dy*.3
     end if


     labarr = new(numbox+2,"string")         ;Prepare actual X-axis labels...
     labarr(0) = ""
     labarr(numbox+1) = ""
     labarr(1:numbox) = yAxis(1:numbox)

; Whether to maximize plot in frame.
     maximize = get_res_value(plotres,"gsnMaximize",False)

     if (plotres) then
;       print("Plot resources detected, accepting")
        fsatts = getvaratts(plotres)
        do tx = 0,dimsizes(fsatts)-1
           if (fsatts(tx).eq."tmYLLabels") then        ;Special section to test for YLLabels
              if (dimsizes(plotres@tmYLLabels).ne.numbox) then
                     print("boxplot_vertical: Fatal:Number of YB Labels does not match number of boxes, exiting")
                     exit
                  else
                     labarr(1:numbox) = plotres@$fsatts(tx)$
                  end if
                  yllab = plotres@tmYLLabels
                  delete(plotres@tmYLLabels)      ;Delete so tmYLLabels is not used when all
           end if                       ;atts are assigned below...
        end do
        delete(fsatts)
     end if

     plot = create "plot" logLinPlotClass wks
        "trXMinF" : min(x)
        "trXMaxF" : 10.^(floattointeger(log10(max(x)))+1)
        "trXLog"  : True
        "trYMinF" : min(yAxis)-0.5
        "trYMaxF" : max(yAxis)+0.5
        "pmTickMarkDisplayMode" : "Always"
        "tmYLMode"   : "Automatic"
        "tmXTOn"  : True
        "tmYROn"  : True
        "tmXTBorderOn" : True
        "tmYRBorderOn" : True
        "pmTitleDisplayMode": "Always"                   ; allow titles
        "tiMainOn"     : True
        "tiMainString"   : ""
     end create
     if (plotres) then
        attsetvalues(plot,plotres)
     end if

     polyres = True                          ;Set up defaults
     polyres@gsLineColor       = "black"                     ; color of lines
     polyres@gsLineThicknessF  = 1.5                          ; thickness of lines
     polyres@gsLineDashPattern = 0
     if (lineres) then
        fsatts = getvaratts(lineres)
        do tx = 0,dimsizes(fsatts)-1
           polyres@$fsatts(tx)$ = lineres@$fsatts(tx)$
        end do
     end if

     if (boxOpts) then
        if (isatt(boxOpts,"boxColors")) then
           boxcolor = boxOpts@boxColors
           if (dimsizes(boxcolor).eq.1.or.dimsizes(boxcolor).ne.numseries) then
              if (dimsizes(boxcolor).ne.1) then
                 print("boxplot_vertical: warning: Number of input colors must either equal 1 or the number of series ("+numseries+").")
                 print(" Using first specified color only:"+boxcolor(0))
              end if
              polyres@gsLineColor = boxcolor(0)
              cflag = 1
           else
              cflag = 2
           end if
        else
           cflag =1
        end if
     else
        cflag = 1
     end if

     dum = new((/numseries,numbox,9/),graphic)
     do ss = 0,numseries-1
       do gg = 0,numbox-1
  ;        print("ss="+ss+" gg="+gg)
          ff = yAxis(gg+1)
  ;        print(gg+" "+ff)
          if (cflag.eq.2) then
              polyres@gsLineColor = boxcolor(ss)
          end if
          xx = (/x(ss,gg,4),x(ss,gg,4)/)
          yy = (/(ff-(boxWidths(gg)/8.)),(ff+(boxWidths(gg)/8.))/)
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,0) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if

          xx = (/x(ss,gg,3),x(ss,gg,4)/)
          yy = (/ff,ff/)
          polyres@gsLineDashPattern = 1
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,1) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if
          polyres@gsLineDashPattern = 0

          xx = (/x(ss,gg,3),x(ss,gg,3)/)
          yy = (/(ff-(boxWidths(gg)/2.)),(ff+(boxWidths(gg)/2.))/)
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,2) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if

          xx = (/x(ss,gg,1),x(ss,gg,3)/)
          yy = (/(ff-(boxWidths(gg)/2.)),(ff-(boxWidths(gg)/2.))/)
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,3) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if

          xx = (/x(ss,gg,2),x(ss,gg,2)/)
          yy = (/(ff-(boxWidths(gg)/2.)),(ff+(boxWidths(gg)/2.))/)
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,4) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if

          xx = (/x(ss,gg,1),x(ss,gg,3)/)
          yy = (/(ff+(boxWidths(gg)/2.)),(ff+(boxWidths(gg)/2.))/)
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,5) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if

          xx = (/x(ss,gg,1),x(ss,gg,1)/)
          yy = (/(ff-(boxWidths(gg)/2.)),(ff+(boxWidths(gg)/2.))/)
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,6) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if

          xx = (/x(ss,gg,0),x(ss,gg,1)/)
          yy = (/ff,ff/)
          polyres@gsLineDashPattern = 1
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,7) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if
          polyres@gsLineDashPattern = 0

          xx = (/x(ss,gg,0),x(ss,gg,0)/)
          yy = (/(ff-(boxWidths(gg)/8.)),(ff+(boxWidths(gg)/8.))/)
          if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
             dum(ss,gg,8) = gsn_add_polyline(wks,plot,xx,yy,polyres)
          end if
       end do
     end do
     dumname = unique_string("dum")
;     print(dumname)
     plot@$dumname$ = dum
     if (isvar("yllab").eq."True") then
        plotres@tmYLLabels = yllab      ;reassign YLLabels for possible future plots
     end if
     if (maximize) then
        mres = True
        mres@gsnDraw  = False
        mres@gsnFrame = False
        maximize_output(wks,mres)
     end if
     print("returning...")
     return(plot)
end
; ---------------------------------
