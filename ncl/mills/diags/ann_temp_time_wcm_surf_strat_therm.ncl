; Prior to running, do:
;   ncrcat -v TS <h0.files*.nc> <TS_timeseries.nc>
; and have that as your addfile input in ncl.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

host=getenv("HOSTNAME")

if (host.eq."cmip5") then
  in_dir = "/data/waccm/"
else
  in_dir = "/data3/mmills/"
end if

minyear = 1960
maxyear=2099

run_name1="b40.1955-2005.2deg.wcm.avg"
run_name2="b40.rcp8_5.2deg.wcm.003"

filename1=in_dir+run_name1+"/proc/b40.rcp4_5.2deg.wcm.ave.cam2.h0zm.annual.TS.T.1960-2099.nc"
filename2=in_dir+run_name2+"/proc/b40.rcp8_5.2deg.wcm.003.cam2.h0zm.annual.TS.T.2005-2099.nc"

file1=addfile(filename1,"r")
file2=addfile(filename2,"r")

norm_plot = True
norm_plot_yrs_abs = (/1961,1990/)

nyears = maxyear - minyear + 1

gw  =file1->gw
lev =file1->lev
time=file1->time
date=file1->date
year= (date/10000)+0.5
ntimes=dimsizes(date)
y2005=ind(year.eq.2005.5)

TS1 =file1->TS
T1  =file1->T
TS2 =file2->TS
T2  =file2->T

Z3  =file2->Z3

T1_zm =dim_avg_wgt_Wrap(T1,gw,0)
TS1_zm=dim_avg_wgt_Wrap(TS1,gw,0)
T2_zm =T1_zm
TS2_zm=TS1_zm
T2_zm(y2005:,:) =dim_avg_wgt_Wrap(T2,gw,0)
TS2_zm(y2005:)=dim_avg_wgt_Wrap(TS2,gw,0)

Z3_zm =dim_avg_wgt_Wrap(Z3(0,:,:),gw,0)
Z3_zm =Z3_zm / 1000.0 ; convert m to km
Z3_zm@units = "km"

if (norm_plot) then

  T_zm_mean=dim_avg_Wrap(T1_zm(lev|:,time|norm_plot_yrs_abs(0)-minyear:norm_plot_yrs_abs(1)-minyear))
  TS_zm_mean=dim_avg_Wrap(TS1_zm(norm_plot_yrs_abs(0)-minyear:norm_plot_yrs_abs(1)-minyear))

  T1_zm=T1_zm-conform(T1_zm,T_zm_mean,1)
  TS1_zm=TS1_zm-TS_zm_mean
  T2_zm=T2_zm-conform(T2_zm,T_zm_mean,1)
  TS2_zm=TS2_zm-TS_zm_mean
  
end if

T1_smooth=runave(T1_zm(lev|:,time|:),11,0)
T2_smooth=runave(T2_zm(lev|:,time|:),11,0)

filename=  "$IMG_PATH/ann_temp_time_wcm_surf_strat_therm"
print(filename)
wks   = gsn_open_wks ("pdf", filename)

plot1 = new(3,graphic)

fheight = 0.05

pan = True
pan@gsnMaximize = True
;pan@gsnPaperMargin = 0.0          ; expand plot size to maximum
pan@gsnPaperOrientation = "landscape"

res = True
res@gsnDraw = False
res@gsnFrame = False

res@trXMaxF = 2099
res@gsnYRefLine           = 0.              ; reference line   
res@gsnYRefLineDashPattern = 1.
res@xyLineColors   = (/"deepskyblue","palevioletred","Blue","Red"/)
res@xyDashPatterns = (/0,0,0,0/)
res@xyMonoLineThickness = False
res@xyLineThicknesses = (/4,4,2,2/)
 

res@vpWidthF        = 2.9
res@tiMainFontHeightF = fheight
res@tiXAxisString = " "
res@tiYAxisFontHeightF = fheight
res@tiYAxisString = "Temperature anomaly (K)"
res@tiYAxisFontHeightF = fheight

 txres               = True                            
 txres@txFontColor = "Blue"             ; Set the font color
 txres2               = True                            
 txres2@txFontColor = "Red"             ; Set the font color

lev1=3
data=conform_dims((/4,ntimes/),T1_zm(lev|lev1,time|:),1)
data!1="time" 
data&time=time
data(3,:)=T2_zm(lev|lev1,time|:)
data(0,:)=runave(data(2,:),11,0)
data(1,:)=runave(data(3,:),11,0)

;label = sprintf("%3.2g",lev(3))+" hPa"
label = sprintf("%4.3g",Z3_zm(lev1))+" km"
plot1(0) = gsn_csm_xy (wks,year,data,res)
;*****************************************************************
; get some info from the plot object to use in constructing an
; extra plot string
;*****************************************************************
 getvalues plot1(0)                     
   "tmYLLabelFontHeightF"   : fheight
   "tmXTValues"             : tmXTValues
   "tmYLValues"             : tmYLValues
 end getvalues
 nTm  = dimsizes(tmXTValues)               ; number of major tick marks
;gsn_text(wks,plot1(0),label,0.75*tmXTValues(nTm-1), \
;                                    0.35*tmYLValues(nTm-1), False)
text=gsn_add_text(wks,plot1(0),label,2080.0, 10.0, False)
text2=gsn_add_text(wks,plot1(0),"RCP4.5",1975.0,-45.0, txres)
text2=gsn_add_text(wks,plot1(0),"RCP8.5",1975.0,-70.0, txres2)

;label = sprintf("%3.2g",lev({0.95}))+" hPa"
label = sprintf("%3.2g",Z3_zm({0.95}))+" km"
data(2,:)=T1_zm({lev|0.95},time|:)
data(3,:)=T1_zm({lev|0.95},time|:)
data(3,:)=T2_zm({lev|0.95},time|:)
data(0,:)=runave(data(2,:),11,0)
data(1,:)=runave(data(3,:),11,0)
plot1(1) = gsn_csm_xy (wks,year,data,res)
text=gsn_add_text(wks,plot1(1),label,2080.0, -2.0, False)

label = "Surface"
data(2,:)=TS1_zm
data(3,:)=TS1_zm
data(3,:)=TS2_zm
data(0,:)=runave(data(2,:),11,0)
data(1,:)=runave(data(3,:),11,0)
plot1(2) = gsn_csm_xy (wks,year,data,res)
text=gsn_add_text(wks,plot1(2),label,2080.0, 0.7, False)

;pan@gsnFrame       = False
gsn_panel(wks,plot1,(/3,1/),pan)

gsn_define_colormap(wks,"BlWhRe")  ; choose a colormap

res2                     = True          ; plot mods desired
res2@gsnDraw = False
res2@gsnFrame = False

res2@cnFillOn            = True          ; turn on color
res2@gsnSpreadColors     = True          ; use full range of colormap

res2@vpWidthF            = 1.0           ; change aspect ratio of plot
;res2@vpHeightF           = 0.3
res2@sfXArray=year
res2@sfYArray=Z3_zm

res2@cnLevelSelectionMode = "ExplicitLevels"
res2@cnLevels=(/-70,-50,-30,-20,-15,-10,-5,-2,-0.5,0.5,2,5,10,15,20,30,50,70/)
res2@lbAutoManage=False
res2@lbLabelFontHeightF=0.018
res2@trXMinF=1965
res2@trXMaxF=2095
res2@tiYAxisString="Approximate Altitude (km)"

res2@tiMainString="RCP4.5"

plot = gsn_csm_contour(wks,T1_smooth,res2)       ; contour the variable
gsn_panel(wks,plot,(/1,1/),pan)

res2@tiMainString="RCP8.5"
plot2 = gsn_csm_contour(wks,T2_smooth,res2)       ; contour the variable
gsn_panel(wks,plot2,(/1,1/),pan)

end
