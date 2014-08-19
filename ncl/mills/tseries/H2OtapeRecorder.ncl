
rootpath=(/ "/glade/scratch/mmills/archive/", \
            "/glade/scratch/mmills/archive/", \
            "/glade/p/cesm/wawg/runs/cesm1_2/" /)

caseid = (/ "b.e12.B1850C5CN.f19_g16.W5ocn.003", \
            "b.e12.B1850C5CN.f19_g16.W5ocn.L32.005", \
            "b.e12.B1850W5CN.f19_g16.002" /)
            
nCases = dimsizes(caseid)

wks = gsn_open_wks ("pdf","~/images/H2OtapeRecorder")

do i = 0, nCases-1
;  filename=rootpath(i)+caseid(i)+"/atm/proc/h0gm/"+caseid(i)+".cam.h0.10N-10S.nc"
  filename=rootpath(i)+caseid(i)+"/atm/proc/h0gm/"+caseid(i)+".cam.h0.tropicsavg.nc"
  print(i+" opening "+filename)
  a=addfile(filename,"r")
  h2o_a=a->Q
  printVarSummary(h2o_a)
  nTimes=dimsizes(h2o_a&time)
  lev=h2o_a&lev
  lev@units="hPa"
  nlev=dimsizes(lev)
  nYears=nTimes/12
  h2o_4d=reshape(h2o_a,(/nYears,12,nlev/))
  h2o_1y=dim_avg_n(h2o_4d,0)
  h2o_1y!0="time"
  h2o_1y!1="lev"
  h2o_2y=new((/nlev,24/), typeof(h2o_1y))
  h2o_2y(:,0:11)=h2o_1y(lev|:,time|:)
  h2o_2y(:,12:23)=h2o_1y(lev|:,time|:)
  h2o_2y=h2o_2y*1e6*28.964/18.015280 ; kg H2O/kg air -> ppmv
  h2o_2y!0="lev"
  h2o_2y&lev=lev
;  print(h2o_2y({100},:))
  
  
  gsn_define_colormap(wks,"rainbow")                ; choose colormap
  res                      = True                  ; plot mods desired
  res@cnFillOn             = True                  ; turn on color
  res@trYMinF           = 20.0
  res@trYMaxF           = 100.0
  res@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  res@cnMinLevelValF       = 2.0                  ; set min contour level
  res@cnMaxLevelValF       = 4.2                  ; set max contour level
  res@cnLevelSpacingF      = 0.1                  ; set contour spacing
  res@tiYAxisString     ="Pressure (hPa)"
  
  plot = gsn_csm_pres_hgt(wks, h2o_2y, res)
 
  delete(plot)
  delete(h2o_a)
  delete(h2o_4d)
  delete(h2o_1y)
  delete(h2o_2y)
  delete(a)
  delete(filename)
  delete(lev)
end do            
