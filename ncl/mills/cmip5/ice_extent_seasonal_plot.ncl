load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

basepath="/data/waccm/"
cases20c=(/"b40.1955-2005.2deg.wcm.002","b40.1955-2005.2deg.wcm.003","b40.1955-2005.2deg.wcm.004"/)
files=basepath+cases20c+"/proc/iextent_reg_avg_timeseries_sh.dat"
ncol=7
ignorelines=1
iextent_table1=readAsciiTable(files(0),ncol,"float",ignorelines)
;print(dimsizes(iextent_table1))

iextent_table2=readAsciiTable(files(1),ncol,"float",ignorelines)
;print(dimsizes(iextent_table2))

iextent_table12=array_append_record(iextent_table1, iextent_table2, 0)
;print(dimsizes(iextent_table12))

iextent_table3=readAsciiTable(files(2),ncol,"float",ignorelines)
;print(dimsizes(iextent_table3))

iextent_table=array_append_record(iextent_table12, iextent_table3, 0)
;print(dimsizes(iextent_table))

year=iextent_table(:,0)

indrange=ind((year.gt.1979).and.(year.le.2006))

time=year(indrange)
iextentSH=iextent_table(indrange,6)

;print(dimsizes(time))
print(time(0:24)+" "+iextentSH(0:24))

iextentSHmn=reshape(iextentSH,(/dimsizes(iextentSH)/12,12/))
;print(iextentSH(0:11)+" -> "+iextentSHmn(0,0:11))
;print(dimsizes(iextentSHmn))
iextentSHmn!1="month"
iextentSHmn!0="year"
iextentSH20Cavg=dim_avg_Wrap(iextentSHmn(month|:,year|:))
iextentSH20Cstd=dim_stddev_Wrap(iextentSHmn(month|:,year|:))

delete([/iextent_table1,iextent_table2,iextent_table12,iextent_table3,iextent_table,time,year,iextentSH,iextentSHmn,year/])

print(iextentSH20Cavg+"  "+iextentSH20Cstd)

case1850="b40.1850.track1.2deg.wcm.007"
file1850=basepath+case1850+"/proc/iextent_reg_avg_timeseries_sh.dat"
;print((/file1850/))
iextent_table1850=readAsciiTable(file1850,ncol,"float",ignorelines)

idims=dimsizes(iextent_table1850)
;print("iextent_table1850:"+idims)
;print(iextent_table1850(0,:))

year=iextent_table1850(:,0)
indrange=ind((year.ge.296-27*3).and.(year.lt.296))
time1850=year(indrange)
nmonths=dimsizes(time1850)
;print("nmonths:"+nmonths)

iextentSH1850=iextent_table1850(indrange,6)
;print("iextentSH1850:"+dimsizes(iextentSH1850))
print(time1850(0:24)+"   "+iextentSH1850(0:24))

iextentSH1850mn=reshape(iextentSH1850,(/dimsizes(iextentSH1850)/12,12/))
;print(iextentSH1850(0:11)+" -> "+iextentSH1850mn(0,0:11))
iextentSH1850mn!1="month"
iextentSH1850mn!0="year"
iextentSH1850avg=dim_avg_Wrap(iextentSH1850mn(month|:,year|:))
iextentSH1850std=dim_stddev_Wrap(iextentSH1850mn(month|:,year|:))

print(" ")
print(iextentSH1850avg+"  "+iextentSH1850std)

; CCSM4 1979-2005 ice extent (sq m):
ccsm4_20c_mean=(/1.57043e+13,  1.04253e+13,  9.90689e+12,  1.39798e+13,  1.74110e+13, 1.98727e+13,  \
                 2.17365e+13, 2.29591e+13,  2.33879e+13,  2.29711e+13,  2.17322e+13, 1.93954e+13 /)
ccsm4_20c_std =(/5.70452e+11,  5.22612e+11,  5.54890e+11,  4.87493e+11,  4.70718e+11, 4.64405e+11,  \
                 4.64461e+11, 4.45525e+11,  4.68508e+11,  4.38912e+11,  4.62639e+11, 4.94050e+11 /)
ccsm4_20c_mean=ccsm4_20c_mean/1e12 ; (m2 -> 1e6 km2)
ccsm4_20c_std=ccsm4_20c_std/1e12 ; (m2 -> 1e6 km2)

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


data = new((/6,12/),"float")
data(0,:)=(/iextentSH1850avg/)
data(1,:)=(/iextentSH20Cavg/)
data(2,:)=(/ccsm4_1850_mean/)
data(3,:)=(/ccsm4_20c_mean/)
data(4,:)=(/obs_mean/)
data(5,:)=(/hindcast_mean/)

time = ispan(0,11,1)

  wks = gsn_open_wks("pdf","$IMG_PATH/ice_extent_seasonal")                 ; create plot
  
  res                    = True              ; plot mods desired
  res@xyDashPatterns = (/0,1,0,1,0,1/)
  res@xyLineColors = (/"red", "red", "blue", "blue", "black", "black"/)
  res@xyMonoLineColor = False
  res@xyLineThicknessF = 6
  res@trXMaxF = 11
  top_plot = gsn_csm_xy (wks,time,data,res)       ; create line plot
  
delete(wks)
