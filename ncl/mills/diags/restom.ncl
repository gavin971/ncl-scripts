load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

host=getenv("HOSTNAME")

run_name = "f.e11.FW.f19_f19.001"
rootpath="/Volumes/MillsPassport/cesm/"
year0 = 1
nyears=9-year0+1

in_dir =rootpath+run_name+"/h0/"
save_dir =rootpath+run_name+"/analysis/"
savefile=save_dir+run_name+"_RESTOM.nc"
  
image_path=getenv("IMAGE_PATH")
pdffile=image_path+"/"+run_name+"_RESTOM"
print("Creating PDF file "+pdffile+".pdf")
wks   = gsn_open_wks ("pdf", pdffile)

mon_wgt = (/0.08493150770664215,0.07671232521533966,0.08493150770664215,\
            0.08219178020954132,0.08493150770664215,0.08219178020954132,\
            0.08493150770664215,0.08493150770664215,0.08219178020954132,\
            0.08493150770664215,0.08219178020954132,0.08493150770664215/)

infiles = systemfunc("ls -1 "+in_dir+run_name+"*h0.*.nc")
nyears=dimsizes(infiles)/12

if (isfilepresent(savefile)) then

  print("reading "+savefile)
  b = addfile(savefile,"r")
  year_a=b->year
  restom_ann_w_a=b->restom_ann_w
  ts_ann_w_a=b->ts_ann_w
  ts_land_ann_w_a=b->ts_land_ann_w
  sst_ann_w_a=b->sst_ann_w

  nyears_a = dimsizes(year_a)
else
  nyears_a=0
end if

if (nyears_a.ge.nyears) then ; use only savefile if it has all history years

  year=year_a
  restom_ann_w=restom_ann_w_a
  ts_ann_w=ts_ann_w_a
  ts_land_ann_w=ts_land_ann_w_a
  sst_ann_w=sst_ann_w_a  
  
else

  i=nyears_a*12
  fin=addfile(infiles(i),"r")
  print((/infiles(i)/))
;  print(fin->date)
  fsnt = fin->FSNT
  flnt = fin->FLNT
  ts = fin->TS
  landfrac = fin->LANDFRAC
  year_in = fin->date ; get record coordinate
  year_in = year_in/10000
;  print(year_in)

  do i=nyears_a*12+1,dimsizes(infiles)-1,1
;    print(i+" "+(/infiles(i)/))
    fin=addfile(infiles(i),"r")
;    print(fin->date)

    temp=fsnt
    delete(fsnt)
    fsnt = array_append_record(temp,fin->FSNT,0)
    delete(temp)

    temp=flnt
    delete(flnt)
    flnt = array_append_record(temp,fin->FLNT,0)
    delete(temp)

    temp=ts
    delete(ts)
    ts = array_append_record(temp,fin->TS,0)
    delete(temp)

    temp=landfrac
    delete(landfrac)
    landfrac = array_append_record(temp,fin->LANDFRAC,0)
    delete(temp)

    temp=year_in
    delete(year_in)
;    print(i+" "+temp+" "+fin->date+" "+(fin->date)/10000)
    year_in = array_append_record(temp,(fin->date)/10000,0)
    delete(temp)

  end do
  ;    fid=addfile(in_file,"r")
  ;    data = fid->$variable$
  ;inptr = addfiles(infiles,"r")
  ;fsnt    = addfiles_GetVar(inptr,infiles,"FSNT")
  ;flnt    = addfiles_GetVar(inptr,infiles,"FLNT")

  ;ts    = addfiles_GetVar(inptr,infiles,"TS")
  ;landfrac = addfiles_GetVar(inptr,infiles,"LANDFRAC")

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

  year = new(nyears,integer)

  if (nyears_a.gt.0) then
      restom_ann_w(0:nyears_a-1) = restom_ann_w_a
      ts_ann_w(0:nyears_a-1) = ts_ann_w_a
      sst_ann_w(0:nyears_a-1) = sst_ann_w_a
      ts_land_ann_w(0:nyears_a-1) = ts_land_ann_w_a
      year(0:nyears_a-1) = year_a
  end if
  
  do iyr=nyears_a,nyears-1
    print((/iyr/))
    iyr2=iyr-nyears_a
    year(iyr)=year_in(iyr2*12)
    do imonth=0,11
      restom_ann_w(iyr) = mon_wgt(imonth)*dim_avg(restom_w(iyr2*12+imonth))+restom_ann_w(iyr)
      ts_ann_w(iyr) = mon_wgt(imonth)*dim_avg(ts_w(iyr2*12+imonth))+ts_ann_w(iyr)
      sst_ann_w(iyr) = mon_wgt(imonth)*dim_avg(sst_w(iyr2*12+imonth))+sst_ann_w(iyr)
      ts_land_ann_w(iyr) = mon_wgt(imonth)*dim_avg(ts_land_w(iyr2*12+imonth))+ts_land_ann_w(iyr)
    end do
  end do
  
  restom_ann_w!0="year"
  restom_ann_w&year=year
  
  ts_ann_w!0="year"
  ts_ann_w&year=year
  
  sst_ann_w!0="year"
  sst_ann_w&year=year
  
  ts_land_ann_w!0="year"
  ts_land_ann_w&year=year

;  year = ispan (year0,nyears+year0-1,1)  

  a = systemfunc("mkdir -p "+save_dir)
  ; Remove any preexisting file.
  print("  Deleting " + savefile)
  system("/bin/rm " + savefile)
  
  print("  Creating " + savefile)
  b = addfile(savefile,"c")
  b->year=year
  b->restom_ann_w=restom_ann_w
  b->ts_ann_w=ts_ann_w
  b->ts_land_ann_w=ts_land_ann_w
  b->sst_ann_w=sst_ann_w

end if

;print(restom_ann_w)

;print("Average = "+dim_avg(restom_ann_w))


; Averages for annual means.

  plot = new(4, graphic)

resP                       = True
resP@gsnMaximize           = True            ; expand plot size to maximum
resP@gsnPaperOrientation   = "landscape"
resP@txString   = run_name


res = True
res@gsnDraw               = False
res@gsnFrame           = False   ; Don't advance the frame

;res@trYMinF = -4
;res@trYMaxF = 4
res@gsnYRefLine           = 0.              ; reference line   
res@gsnXYBarChart         = True            ; create bar chart
res@gsnAboveYRefLineColor = "red"           ; above ref line fill red
res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
res@tiYAxisString = "RESTOM (W/m2)" ; y-axis label      
res@tiMainString  = run_name ; title


  res1 = True   ; These four lines are redundant.
  res2 = True
  res3 = True
  res4 = True

  res1 = res
  res2 = res
  res3 = res
  res4 = res


print(year+" "+restom_ann_w+" "+ts_ann_w+" "+ts_land_ann_w+" "+sst_ann_w)

res1@tiMainString  = "RESTOM" ; title
res1@gsnYRefLine           = 0.              ; reference line   
;res1@trYMinF = -4
res1@trYMaxF = 6
res1@vpXF              = 0.1    ; X,Y position of first plot
res1@vpYF              = 0.91
print("year:"+dimsizes(year))
print("restom_ann_w:"+dimsizes(restom_ann_w))
print(min(restom_ann_w)+" "+max(restom_ann_w))
plot(0) = gsn_csm_xy (wks,year,restom_ann_w,res1)

res2@tiMainString  = "TS" ; title
res2@tiYAxisString = "TS K" ; y-axis label      
res2@vpXF              = 0.6    ; X,Y position of second plot
res2@vpYF              = 0.91
;res2@trYMinF = 286.0
res2@trYMinF = 286.0
res2@gsnYRefLine           =  res2@trYMinF             ; reference line   
;res2@trYMaxF = 289
res2@trYMaxF = 290.0  ;caf
print("ts_ann_w:"+dimsizes(ts_ann_w))
print(min(ts_ann_w)+" "+max(ts_ann_w))
print(res2)
plot(1) = gsn_csm_xy (wks,year,ts_ann_w,res2)

res3@tiMainString  = "TS_LAND" ; title
res3@vpXF              = 0.1    ; X,Y position of third plot
res3@vpYF              = 0.45
;res3@trYMinF = 280.0
res3@trYMinF = 280.0
res3@gsnYRefLine           = res3@trYMinF              ; reference line   
;res3@trYMaxF = 283
res3@trYMaxF = 284.0  ;caf
res3@tiYAxisString = "TS_LAND K" ; y-axis label      
print("ts_land_ann_w:"+dimsizes(ts_land_ann_w))
print(min(ts_land_ann_w)+" "+max(ts_land_ann_w))
plot(2) = gsn_csm_xy (wks,year,ts_land_ann_w,res3)

res4@tiMainString  = "SST" ; title
res4@vpXF              = 0.6    ; X,Y position of fourth plot
res4@vpYF              = 0.45
res4@tiYAxisString = "SST C" ; y-axis label      
res4@trYMinF = 16
res4@trYMaxF = 19
print("sst_ann_w:"+dimsizes(sst_ann_w))
print(min(sst_ann_w)+" "+max(sst_ann_w))
plot(3) = gsn_csm_xy (wks,year,sst_ann_w,res4)

gsn_panel(wks, plot, (/ 2, 2 /), resP)


delete(wks)

print("gv "+pdffile+".pdf")
print("mutt -a "+pdffile+".pdf mikejmills@mac.com < /dev/null")

;mutter=systemfunc("mutt -a " + pdffile + ".pdf mikejmills@mac.com < /dev/null")
;print(mutter)

