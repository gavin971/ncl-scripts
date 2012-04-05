load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


;begin

;year0 = 96
;run_name = "b40.1850.track1.2deg.wcm.007"

year0 = 1
year1 = 8
nyears=year1-year0+1

run_name=new(2,string)
rootpath=new(2,string)

run_name(0) = "f.e11.FW5.f19_f19.001"
rootpath(0) = "/Volumes/MillsPassport/cesm/"

run_name(1) = "f.e11.FW.f19_f19.001"
rootpath(1) = "/Volumes/MillsPassport/cesm/"

in_dir   = rootpath+run_name+"/h0/"
save_dir = rootpath+run_name+"/analysis/"
savefile = save_dir+run_name+"_RESTOM.nc"

image_path=getenv("IMAGE_PATH")
pdffile=image_path+"/"+run_name(0)+"-"+run_name(1)+"_RESTOM"
print("Creating PDF file "+pdffile+".pdf")
wks   = gsn_open_wks ("pdf", pdffile)

mon_wgt = (/0.08493150770664215,0.07671232521533966,0.08493150770664215,\
            0.08219178020954132,0.08493150770664215,0.08219178020954132,\
            0.08493150770664215,0.08493150770664215,0.08219178020954132,\
            0.08493150770664215,0.08219178020954132,0.08493150770664215/)

do irun=0,1
;  print(irun)
  if (irun.eq.1) then
    delete(infiles)
    delete(b)
    delete(year_a)
    delete(restom_ann_w_a)
    delete(ts_ann_w_a)
    delete(ts_land_ann_w_a)
    delete(sst_ann_w_a)
    delete(nyears_a)
  end if
  infiles = systemfunc("ls -1 "+in_dir(irun)+run_name(irun)+"*h0.*nc")
  nyears=dimsizes(infiles)/12

  if (isfilepresent(savefile(irun))) then
    
    b = addfile(savefile(irun),"r")
    year_a=b->year
    restom_ann_w_a=b->restom_ann_w
    ts_ann_w_a=b->ts_ann_w
    ts_land_ann_w_a=b->ts_land_ann_w
    sst_ann_w_a=b->sst_ann_w

    nyears_a = dimsizes(year_a)
;    print(year_a)
  else
    print ("file missing: "+savefile(irun))
    return
  end if

  if (nyears_a.ge.nyears) then

    if (irun.eq.0) then
      year=year_a
      restom_ann_w=restom_ann_w_a
      ts_ann_w=ts_ann_w_a
      ts_land_ann_w=ts_land_ann_w_a
      sst_ann_w=sst_ann_w_a  
    else ; irun.eq.1
      year2=year_a
      restom_ann_w2=restom_ann_w_a
      ts_ann_w2=ts_ann_w_a
      ts_land_ann_w2=ts_land_ann_w_a
      sst_ann_w2=sst_ann_w_a  
    end if ; irun.eq.0
  else;(nyears_a.lt.nyears)

    print ("missing years in "+savefile(irun))
    print ("nyears="+nyears+", nyears_a="+nyears_a)
    return

  end if ; (nyears_a.ge.nyears)
end do ; irun=0,1

;print(restom_ann_w)

;print("Average = "+dim_avg(restom_ann_w))


; Averages for annual means.


res = True
res@gsnFrame           = False   ; Don't advance the frame
res@vpHeightF          = 0.32    ; Width and height of
res@vpWidthF           = 0.32    ; each plot in NDC space.

;res@trYMinF = -4
;res@trYMaxF = 4
res@gsnYRefLine           = 0.              ; reference line   
res@gsnXYBarChart         = False           ; create line chart
;res@gsnAboveYRefLineColor = "red"           ; above ref line fill red
;res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
res@xyLineColor = "red"
res@xyMonoLineColor = True
res@xyLineThicknessF = 2.0
res@xyMonoLineThickness = True
res@tiYAxisString = "RESTOM (W/m2)" ; y-axis label      
res@tiMainString  = run_name(0)+"-"+run_name(1) ; title
res@trXMinF = year0
res@trXMaxF = year0+nyears

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
gsn_text_ndc(wks,run_name(0)+"-"+run_name(1),0.5,.98,txres)

print(year+" "+restom_ann_w+" "+ts_ann_w+" "+ts_land_ann_w+" "+sst_ann_w)

res1@tiMainString  = "RESTOM" ; title
res1@gsnYRefLine           = 0.              ; reference line   
res1@trYMinF = -0.5
res1@trYMaxF = 2
res1@vpXF              = 0.1    ; X,Y position of first plot
res1@vpYF              = 0.91
res1black=res1
res1black@xyLineColor = "black"
print("year:"+dimsizes(year))
print("restom_ann_w:"+dimsizes(restom_ann_w))
print(min(restom_ann_w)+" "+max(restom_ann_w))
plot = gsn_csm_xy (wks,year,restom_ann_w,res1)
plot = gsn_csm_xy (wks,year2,restom_ann_w2,res1black)

res2@tiMainString  = "TS" ; title
res2@tiYAxisString = "TS K" ; y-axis label      
res2@vpXF              = 0.6    ; X,Y position of second plot
res2@vpYF              = 0.91
;res2@trYMinF = 286.0
res2@trYMinF = 288.5
res2@gsnYRefLine           =  res2@trYMinF             ; reference line   
;res2@trYMaxF = 289
res2@trYMaxF = 290.0  ;caf
res2black=res2
res2black@xyLineColor = "black"
print("ts_ann_w:"+dimsizes(ts_ann_w))
print(min(ts_ann_w)+" "+max(ts_ann_w))
;print(res2)
plot = gsn_csm_xy (wks,year,ts_ann_w,res2)
plot = gsn_csm_xy (wks,year2,ts_ann_w2,res2black)

res3@tiMainString  = "TS_LAND" ; title
res3@vpXF              = 0.1    ; X,Y position of third plot
res3@vpYF              = 0.45
;res3@trYMinF = 280.0
res3@trYMinF = 282.0
res3@gsnYRefLine           = res3@trYMinF              ; reference line   
;res3@trYMaxF = 283
res3@trYMaxF = 284.5  ;caf
res3@tiYAxisString = "TS_LAND K" ; y-axis label      
res3black=res3
res3black@xyLineColor = "black"
print("ts_land_ann_w:"+dimsizes(ts_land_ann_w))
print(min(ts_land_ann_w)+" "+max(ts_land_ann_w))
plot = gsn_csm_xy (wks,year,ts_land_ann_w,res3)
plot = gsn_csm_xy (wks,year2,ts_land_ann_w2,res3black)

res4@tiMainString  = "SST" ; title
res4@vpXF              = 0.6    ; X,Y position of fourth plot
res4@vpYF              = 0.45
res4@tiYAxisString = "SST C" ; y-axis label      
res4@trYMinF = 18.0
res4@trYMaxF = 19.5
res4black=res4
res4black@xyLineColor = "black"
print("sst_ann_w:"+dimsizes(sst_ann_w))
print(min(sst_ann_w)+" "+max(sst_ann_w))
plot = gsn_csm_xy (wks,year,sst_ann_w,res4)
plot = gsn_csm_xy (wks,year2,sst_ann_w2,res4black)

print("gv "+pdffile+".pdf")
print("mutt -s RESTOM -a "+pdffile+".pdf mikejmills@mac.com < /dev/null")

frame(wks)
;opentest = systemfunc("open "+pdffile+".pdf")
;mutter=systemfunc("mutt -s RESTOM -a "+pdffile+".pdf mikejmills@mac.com < /dev/null")
;print(mutter)

;end

