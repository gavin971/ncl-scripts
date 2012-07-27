load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

year_bounds=(/(/1960, 1979/),(/1979,2005/)/)

; WACCM4 1979-2005 ice extent
basepath="/data/waccm/"
cases20c=(/"b40.1955-2005.2deg.wcm.002","b40.1955-2005.2deg.wcm.003","b40.1955-2005.2deg.wcm.004"/)
ncases=dimsizes(cases20c)
delete(files)
files=basepath+cases20c+"/proc/iextent_reg_avg_timeseries_sh.dat"
ncol=7
ignorelines=1
print("reading "+files(0))
iextent_table=readAsciiTable(files(0),ncol,"float",ignorelines)

do i = 1, ncases-1
  print("reading "+files(i))
  iextent_table2=readAsciiTable(files(i),ncol,"float",ignorelines)
  iextent_table12=array_append_record(iextent_table, iextent_table2, 0)
  delete(iextent_table)
  iextent_table=iextent_table12
  delete(iextent_table12)
  delete(iextent_table2)
end do

year=iextent_table(:,0)

iextentSH20Cavg=new((/2,12/),"float")
iextentSH20Cstd=new((/2,12/),"float")
do i = 0, 1
  indrange=ind((year.gt.year_bounds(i,0)).and.(year.le.year_bounds(i,1)+1))

  time=year(indrange)
  iextentSH=iextent_table(indrange,6)

  ;print(dimsizes(time))
  print("time(0:24) iextentSH(0:24):")
  print(time(0:24)+" "+iextentSH(0:24))

  iextentSHmn=reshape(iextentSH,(/dimsizes(iextentSH)/12,12/))
  ;print(iextentSH(0:11)+" -> "+iextentSHmn(0,0:11))
  ;print(dimsizes(iextentSHmn))
  iextentSHmn!1="month"
  iextentSHmn!0="year"
  iextentSH20Cavg(i,:)=dim_avg_Wrap(iextentSHmn(month|:,year|:))
  iextentSH20Cstd(i,:)=dim_stddev_Wrap(iextentSHmn(month|:,year|:))
  delete(indrange)
  delete(time)
  delete(iextentSH)
  delete(iextentSHmn)
end do

print("iextentSH20Cavg  iextentSH20Cstd:")
print(iextentSH20Cavg(0,:)+"  "+iextentSH20Cstd(0,:)+"  "+iextentSH20Cavg(1,:)+"  "+iextentSH20Cstd(1,:))

; WACCM4 pre-industrial ice extent
case1850="b40.1850.track1.2deg.wcm.007"
file1850=basepath+case1850+"/proc/iextent_reg_avg_timeseries_sh.dat"
;print((/file1850/))
iextent_table1850=readAsciiTable(file1850,ncol,"float",ignorelines)

idims=dimsizes(iextent_table1850)
;print("iextent_table1850:"+idims)
;print(iextent_table1850(0,:))

year1850=iextent_table1850(:,0)
indrange1850=ind((year1850.ge.296-27*3).and.(year1850.lt.296))
time1850=year1850(indrange1850)
nmonths=dimsizes(time1850)
;print("time1850(0::12): ")
;print(time1850(0::12))
;print("nmonths:"+nmonths)

iextentSH1850=iextent_table1850(indrange1850,6)
;print("iextentSH1850:"+dimsizes(iextentSH1850))
print("time1850(0:24)   iextentSH1850(0:24):")
print(time1850(0:24)+"   "+iextentSH1850(0:24))

iextentSH1850mn=reshape(iextentSH1850,(/dimsizes(iextentSH1850)/12,12/))
;print(iextentSH1850(0:11)+" -> "+iextentSH1850mn(0,0:11))
iextentSH1850mn!1="month"
iextentSH1850mn!0="year"
iextentSH1850avg=dim_avg_Wrap(iextentSH1850mn(month|:,year|:))
iextentSH1850std=dim_stddev_Wrap(iextentSH1850mn(month|:,year|:))

time1850mn=reshape(time1850,(/dimsizes(iextentSH1850)/12,12/))
time1850mn!1="month"
time1850mn!0="year"
time1850avg=dim_avg_Wrap(time1850mn(month|:,year|:))

print("time1850avg  iextentSH1850avg  iextentSH1850std:")
print(time1850avg+"  "+iextentSH1850avg+"  "+iextentSH1850std)

; CCSM4 1979-2005 ice extent
basepath="/data/ccsm4/"
delete(cases20c)
cases20c=(/"b40.20th.track1.1deg.005","b40.20th.track1.1deg.006","b40.20th.track1.1deg.007",\
           "b40.20th.track1.1deg.008","b40.20th.track1.1deg.009","b40.20th.track1.1deg.012"/)
ncases=dimsizes(cases20c)
delete(files)
files=basepath+cases20c+"/proc/iextent_reg_avg_timeseries_sh.dat"
ncol=7
ignorelines=1
print("reading "+files(0))
delete(iextent_table)
iextent_table=readAsciiTable(files(0),ncol,"float",ignorelines)

do i = 1, ncases-1
  print("reading "+files(i))
  iextent_table2=readAsciiTable(files(i),ncol,"float",ignorelines)
;  print("iextent_table2(:,0):"+iextent_table2(:,0))
  iextent_table12=array_append_record(iextent_table, iextent_table2, 0)
  delete(iextent_table)
  iextent_table=iextent_table12
  delete(iextent_table12)
  delete(iextent_table2)
end do
;exit

delete(year)
year=iextent_table(:,0)

ccsm4_20C_mean=new((/2,12/),"float")
ccsm4_20C_std=new((/2,12/),"float")
do i = 0, 1
  indrange=ind((year.gt.year_bounds(i,0)).and.(year.le.year_bounds(i,1)+1))

  time=year(indrange)
  iextentSH=iextent_table(indrange,6)

  ;print(dimsizes(time))
  print("time(0:24) iextentSH(0:24):")
  print(time(0:24)+" "+iextentSH(0:24))

  iextentSHmn=reshape(iextentSH,(/dimsizes(iextentSH)/12,12/))
  ;print(iextentSH(0:11)+" -> "+iextentSHmn(0,0:11))
  ;print(dimsizes(iextentSHmn))
  iextentSHmn!1="month"
  iextentSHmn!0="year"
  ccsm4_20C_mean(i,:)=dim_avg_Wrap(iextentSHmn(month|:,year|:))
  ccsm4_20C_std(i,:)=dim_stddev_Wrap(iextentSHmn(month|:,year|:))
  delete(indrange)
  delete(time)
  delete(iextentSH)
  delete(iextentSHmn)
end do

; CCSM4 1979-2005 ice extent (sq m):
;ccsm4_20C_mean=(/1.57043e+13,  1.04253e+13,  9.90689e+12,  1.39798e+13,  1.74110e+13, 1.98727e+13,  \
;                 2.17365e+13, 2.29591e+13,  2.33879e+13,  2.29711e+13,  2.17322e+13, 1.93954e+13 /)
;ccsm4_20C_std =(/5.70452e+11,  5.22612e+11,  5.54890e+11,  4.87493e+11,  4.70718e+11, 4.64405e+11,  \
;                 4.64461e+11, 4.45525e+11,  4.68508e+11,  4.38912e+11,  4.62639e+11, 4.94050e+11 /)
ccsm4_20C_mean=ccsm4_20C_mean/1e12 ; (m2 -> 1e6 km2)
ccsm4_20C_std=ccsm4_20C_std/1e12 ; (m2 -> 1e6 km2)

iextentSH1850avg=iextentSH1850avg/1e12 ; (m2 -> 1e6 km2)
iextentSH1850std=iextentSH1850std/1e12 ; (m2 -> 1e6 km2)
iextentSH20Cavg=iextentSH20Cavg/1e12 ; (m2 -> 1e6 km2)
iextentSH20Cstd=iextentSH20Cstd/1e12 ; (m2 -> 1e6 km2)


; CCSM4  preindustrial ice extent (millions of sq km):
ccsm4_1850_mean= (/17.5771, 12.3631, 11.7589, 15.8442, 19.4190, 22.0232, 23.9711, 25.1854, 25.5522, 25.0638, 23.7774, 21.4125/)
ccsm4_1850_std = (/1.06718, 0.997159, 0.881124, 0.662999, 0.604378, 0.567596, 0.558846, 0.554561, 0.557567, 0.555505, 0.560376, 0.649233/)

; ice-ocean hindcast run (millions of sq km):
hindcast_mean=(/6.79222, 2.05835, 1.88318, 6.04499, 10.8974, 14.0386, 16.4871, 18.1812, 18.9177, 18.6684, 17.0007, 12.9849/)

; observations (millions of sq km)
obs_mean=(/5.08556, 3.02296, 4.36556, 7.33593, 10.7170, 13.8504, 16.3848, 18.0589, 18.7211, 18.2863, 16.2341, 11.1389/)

waccmRatio=iextentSH20Cavg
ccsm4Ratio=ccsm4_20C_mean
do i = 0,1
  waccmRatio(i,:)=iextentSH20Cavg(i,:)/iextentSH1850avg
  ccsm4Ratio(i,:)=ccsm4_20C_mean(i,:)/ccsm4_1850_mean  
end do

print(" ")
print("WACCM ratio, CCSM4 ratio")
print(waccmRatio+"   "+ccsm4Ratio)

data = new((/6,12/),"float")
data(0,:)=(/iextentSH1850avg/)
data(1,:)=(/iextentSH20Cavg(1,:)/)
data(2,:)=(/ccsm4_1850_mean/)
data(3,:)=(/ccsm4_20C_mean(1,:)/)
data(4,:)=(/obs_mean/)
data(5,:)=(/hindcast_mean/)

data2 = new((/4,12/),"float")
data2(0,:)=(/waccmRatio(0,:)/)
data2(1,:)=(/waccmRatio(1,:)/)
data2(2,:)=(/ccsm4Ratio(0,:)/)
data2(3,:)=(/ccsm4Ratio(1,:)/)
data2 = 100.0*(data2-1) ; convert ratio to percent loss

mnmx = new((/8,12/),"float")
mnmx(0,:)=(/iextentSH1850avg+iextentSH1850std/)
mnmx(1,:)=(/iextentSH1850avg-iextentSH1850std/)
mnmx(2,:)=(/iextentSH20Cavg(1,:)+iextentSH20Cstd(1,:)/)
mnmx(3,:)=(/iextentSH20Cavg(1,:)-iextentSH20Cstd(1,:)/)
mnmx(4,:)=(/ccsm4_1850_mean+ccsm4_1850_std/)
mnmx(5,:)=(/ccsm4_1850_mean-ccsm4_1850_std/)
mnmx(6,:)=(/ccsm4_20C_mean(1,:)+ccsm4_20C_std(1,:)/)
mnmx(7,:)=(/ccsm4_20C_mean(1,:)-ccsm4_20C_std(1,:)/)

time = ispan(0,11,1)
imagefile="$IMG_PATH/ice_extent_seasonal_"+year_bounds(1,0)+"-"+year_bounds(1,1)

  wks = gsn_open_wks("eps",imagefile)                 ; create plot
  
  res                    = True              ; plot mods desired
  res@gsnDraw            = False             ; don't draw yet
  res@gsnFrame           = False             ; don't advance frame yet
  res@xyDashPatterns = (/1,0,1,0,0,1/)
  res@xyLineColors = (/"red", "red", "blue", "blue", "black", "black"/)
  res@xyMonoLineColor = False
  res@xyLineThicknessF = 4
  res@trXMaxF = 11
  res@trYMinF = 0
  res@trYMaxF = 27
  res@tmXTMinorOn = False
  res@tmXBLabels =(/"J","F","M","A","M","J","J","A","S","O","N","D"/)
  res@tmXBMode   = "Explicit"         ; Define own tick mark labels.
  res@tmXBValues = ispan(0,11,1)
  
  plot = new(2,graphic)
  
  top_plot = gsn_csm_xy (wks,time,data,res)       ; create line plot
  
; Create a plot with the area between both curves filled in blue.
  res2=res
  delete(res2@xyLineColors)
  res2@gsnXYFillColors = "LightPink"
  res2@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
  plot1  = gsn_csm_xy (wks,time,mnmx(0:1,:),res2)  ; Create filled XY plot.
  plot2  = gsn_csm_xy (wks,time,mnmx(2:3,:),res2)  ; Create filled XY plot.

; Create a plot with the area between both curves filled in pink.
  res2@gsnXYFillColors = "LightBlue"
  res2@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
  plot3  = gsn_csm_xy (wks,time,mnmx(4:5,:),res2)  ; Create another filled XY plot.
  plot(0)  = gsn_csm_xy (wks,time,mnmx(6:7,:),res2)  ; Create another filled XY plot.

;
; Overlay the top and mid plots on the bottom plot.
;
; Don't draw anything yet, because we still need to
; attach a legend and some titles.
;
  overlay(plot(0),plot3)
  overlay(plot(0),plot2)
  overlay(plot(0),plot1)
  overlay(plot(0),top_plot)
  
  delete(res@trYMinF)
  delete(res@trYMaxF)
  print(time)
  print(data2)
  plot(1) = gsn_csm_xy (wks,time,data2,res)       ; create line plot


  resP                            = True
  resP@gsnMaximize                = True
  resP@gsnPanelXWhiteSpacePercent = 3
  resP@gsnPaperOrientation        = "landscape"
  resP@gsnPanelLabelBar           = True
  
  gsn_panel(wks,plot,(/1,2/),resP)
;  frame(wks)
  
  delete(wks)
  
  muttcmd="mutt mikejmills@mac.com -a "+imagefile+".eps -s "+imagefile+" < /dev/null"
  print(muttcmd)
  mutter=systemfunc(muttcmd)
  print(mutter)
  print("gv "+imagefile+".eps")
