load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


begin





year0 = 2
nyears = 10
run_name = "b40.1850.waccm.rhminl917"

;year0 = 1
;nyears = 4
;run_name = "b40.1850.waccm.001.test"

;year0 = 1
;nyears = 5
;run_name = "b40.1850.waccm.924vfall"

;year0 = 1
;nyears = 20
;run_name = "b40.1850.2deg.tr1.waccm.001.rhminl924"

;year0 = 1
;nyears = 2
;run_name = "b40.1850.2deg.tr1.waccm.001.5vfall"

;year0 = 1
;nyears = 20
;run_name = "b40.1850.2deg.tr1.waccm.001"

year0 = 1
nyears = 26
run_name = "b40.1850.waccm.rhminl910"

year0 = 40
nyears = 20
run_name = "b40.1850.track1.009"

year0 = 20
nyears = 30
run_name = "b40.1850.track1.2deg.wcm.003"

year0 = 1974
nyears = 30
run_name = "b40.20th.track1.2deg.wcm.005"

year0 = 80
nyears = 10
run_name = "b40.1850.track1.2deg.wcm.004"

year0 = 80
nyears = 10
run_name = "b40.1850.track1.2deg.wcm.003"

in_dir ="/data1/fischer/"+run_name+"/"


mon_wgt = (/0.08493150770664215,0.07671232521533966,0.08493150770664215,0.08219178020954132,0.08493150770664215,0.08219178020954132,0.08493150770664215,0.08493150770664215,0.08219178020954132,0.08493150770664215,0.08219178020954132,0.08493150770664215/)

wks   = gsn_open_wks ("ps", run_name+"_trend")


;infiles = systemfunc("ls -1 "+in_dir+"h0/"+run_name+"*h0.*nc")
infiles = systemfunc("ls -1 "+in_dir+"h0/"+run_name+"*h0.008*nc")
inptr = addfiles(infiles,"r")
fsnt    = addfiles_GetVar(inptr,infiles,"FSNT")
flnt    = addfiles_GetVar(inptr,infiles,"FLNT")

ts    = addfiles_GetVar(inptr,infiles,"TS")
landfrac = addfiles_GetVar(inptr,infiles,"LANDFRAC")
ts_land = mask(ts,landfrac.gt.0.5,True)
sst = mask(ts,landfrac.le.0.33,True)-273




inptr_gw = addfile(infiles(0),"r")
gw1  = inptr_gw->gw

restom = fsnt
restom = (/fsnt-flnt/)


restom_w = wgt_areaave(restom,gw1,1.0,0) 
restom_ann_w = new(nyears,typeof(restom_w))
restom_ann_w = 0.

ts_w = wgt_areaave(ts,gw1,1.0,0) 
ts_ann_w = new(nyears,typeof(ts_w))
ts_ann_w = 0.

 sst_w = wgt_areaave(sst,gw1,1.0,0) 
 sst_ann_w = new(nyears,typeof(sst_w))
 sst_ann_w = 0.

ts_land_w = wgt_areaave(ts_land,gw1,1.0,0) 
ts_land_ann_w = new(nyears,typeof(ts_land_w))
ts_land_ann_w = 0.

do iyr=0,nyears-1
  print(iyr)

  do imonth=0,11
    restom_ann_w(iyr) = mon_wgt(imonth)*dim_avg(restom_w(iyr*12+imonth))+restom_ann_w(iyr)
    ts_ann_w(iyr) = mon_wgt(imonth)*dim_avg(ts_w(iyr*12+imonth))+ts_ann_w(iyr)
    sst_ann_w(iyr) = mon_wgt(imonth)*dim_avg(sst_w(iyr*12+imonth))+sst_ann_w(iyr)
    ts_land_ann_w(iyr) = mon_wgt(imonth)*dim_avg(ts_land_w(iyr*12+imonth))+ts_land_ann_w(iyr)
  end do
end do


;print(restom_ann_w)

;print("Average = "+dim_avg(restom_ann_w))


; Averages for annual means.


res = True
res@gsnFrame           = False   ; Don't advance the frame
  res@vpHeightF          = 0.32    ; Width and height of
  res@vpWidthF           = 0.32    ; each plot in NDC space.

res@trYMinF = -4
res@trYMaxF = 4
res@gsnYRefLine           = 0.              ; reference line   
res@gsnXYBarChart         = True            ; create bar chart
res@gsnAboveYRefLineColor = "red"           ; above ref line fill red
res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
res@tiYAxisString = "RESTOM (W/m2)" ; y-axis label      
res@tiMainString  = run_name ; title


year = ispan (year0,nyears+year0-1,1)  

  res1 = True   ; These four lines are redundant.
  res2 = True
  res3 = True
  res4 = True

  res1 = res
  res2 = res
  res3 = res
  res4 = res


 txres               = True   ; Set up resource list for text strings
  txres@txFontHeightF = 0.025
  gsn_text_ndc(wks,run_name,0.5,.98,txres)


res1@tiMainString  = "RESTOM" ; title
res1@vpXF              = 0.1    ; X,Y position of first plot
res1@vpYF              = 0.91
plot = gsn_csm_xy (wks,year,restom_ann_w,res1)

res2@tiMainString  = "TS" ; title
res2@tiYAxisString = "TS K" ; y-axis label      
res2@vpXF              = 0.6    ; X,Y position of second plot
res2@vpYF              = 0.91
res2@trYMinF = 287
res2@trYMaxF = 289
res2@trYMaxF = 290  ;caf
plot = gsn_csm_xy (wks,year,ts_ann_w,res2)

res3@tiMainString  = "TS_LAND" ; title
res3@vpXF              = 0.1    ; X,Y position of third plot
res3@vpYF              = 0.45
res3@trYMinF = 281
res3@trYMaxF = 283
res3@trYMaxF = 284  ;caf
res3@tiYAxisString = "TS_LAND K" ; y-axis label      
plot = gsn_csm_xy (wks,year,ts_land_ann_w,res3)

res4@tiMainString  = "SST" ; title
res4@vpXF              = 0.6    ; X,Y position of fourth plot
res4@vpYF              = 0.45
res4@tiYAxisString = "SST C" ; y-axis label      
res4@trYMinF = 17
res4@trYMaxF = 19
plot = gsn_csm_xy (wks,year,sst_ann_w,res4)

frame(wks)
end

