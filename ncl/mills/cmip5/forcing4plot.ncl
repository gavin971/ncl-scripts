load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

rootdir="/data/waccm/forcing/"
flnm26="LBC_1765-2500_1.9x2.5_CMIP5_RCP3PD_za_091214.nc"
flnm45="LBC_1765-2500_1.9x2.5_CMIP5_RCP45_za_c091214.nc"
flnm85="LBC_1765-2500_1.9x2.5_CMIP5_RCP85_za_c091214.nc"

rcp26=addfile(rootdir+flnm26,"r")
rcp45=addfile(rootdir+flnm45,"r")
rcp85=addfile(rootdir+flnm85,"r")
gwfile=addfile("/data/waccm/b40.1955-2005.2deg.wcm.avg/proc/"+\
       "b40.1955-2005.2deg.wcm.002.cam2.h0.gw.20th.nc","r")

date=rcp26->date
time=rcp26->time
lat=rcp26->lat  
gw=gwfile->gw

CFC_26=rcp26->CFCL3_LBC
CFC_26=CFC_26+rcp26->CF2CL2_LBC
CO2_26=rcp26->CO2_LBC
CH4_26=rcp26->CH4_LBC
N2O_26=rcp26->N2O_LBC

CFC_45=rcp45->CFCL3_LBC
CFC_45=CFC_45+rcp45->CF2CL2_LBC
CO2_45=rcp45->CO2_LBC
CH4_45=rcp45->CH4_LBC
N2O_45=rcp45->N2O_LBC

CFC_85=rcp85->CFCL3_LBC
CFC_85=CFC_85+rcp85->CF2CL2_LBC
CO2_85=rcp85->CO2_LBC
CH4_85=rcp85->CH4_LBC
N2O_85=rcp85->N2O_LBC

year_start=1850
year_end=2100
nyears=year_end-year_start+1

cfc=new((/3,nyears/),float)
co2=cfc
ch4=cfc
n2o=cfc
cfc@long_name="CFC11+CFC12"
co2@long_name="CO2"
ch4@long_name="CH4"
n2o@long_name="N2O"

do year=year_start,year_end
  print((/year/))
  ydate1=year*10000+101
  ydate2=ydate1+10000
  yind1=closest_val(ydate1,date)
  yind2=closest_val(ydate2,date)-1
  yind3=year-year_start
  
  cfc(0,yind3)=dim_avg(dim_avg_wgt(CFC_26(yind1:yind2,:),gw,1))
  co2(0,yind3)=dim_avg(dim_avg_wgt(CO2_26(yind1:yind2,:),gw,1))
  ch4(0,yind3)=dim_avg(dim_avg_wgt(CH4_26(yind1:yind2,:),gw,1))
  n2o(0,yind3)=dim_avg(dim_avg_wgt(N2O_26(yind1:yind2,:),gw,1))
  
  cfc(1,yind3)=dim_avg(dim_avg_wgt(CFC_45(yind1:yind2,:),gw,1))
  co2(1,yind3)=dim_avg(dim_avg_wgt(CO2_45(yind1:yind2,:),gw,1))
  ch4(1,yind3)=dim_avg(dim_avg_wgt(CH4_45(yind1:yind2,:),gw,1))
  n2o(1,yind3)=dim_avg(dim_avg_wgt(N2O_45(yind1:yind2,:),gw,1))
  
  cfc(2,yind3)=dim_avg(dim_avg_wgt(CFC_85(yind1:yind2,:),gw,1))
  co2(2,yind3)=dim_avg(dim_avg_wgt(CO2_85(yind1:yind2,:),gw,1))
  ch4(2,yind3)=dim_avg(dim_avg_wgt(CH4_85(yind1:yind2,:),gw,1))
  n2o(2,yind3)=dim_avg(dim_avg_wgt(N2O_85(yind1:yind2,:),gw,1))

end do

co2=co2*1e6 ;ppmv
co2@units="ppmv"
ch4=ch4*1e6 ;ppmv
ch4@units="ppmv"
n2o=n2o*1e9 ;ppbv
n2o@units="ppbv"
cfc=cfc*1e12 ;pptv
cfc@units="pptv"

delete(year)
year=ispan(year_start,year_end,1)
year@long_name="year"

wks = gsn_open_wks("ps", "$IMG_PATH/forcing4plot")

plot = new(4, graphic)

resP=True
  resP@gsnPaperOrientation = "landscape"
  resP@gsnMaximize          = True

res=True
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@gsnMaximize          = True
  res@xyLineColors      = (/"blue","green","red"/)
  res@xyDashPatterns      = (/0,0,0/)
  res@vpHeightF             = 0.7              ; change aspect ratio
  res@vpWidthF              = 1.0
  res@trXMinF = 1850
  res@trXMaxF = 2100

  plot(0) = gsn_csm_xy(wks,year,co2,res)
  plot(1) = gsn_csm_xy(wks,year,ch4,res)
  plot(2) = gsn_csm_xy(wks,year,n2o,res)
  plot(3) = gsn_csm_xy(wks,year,cfc,res)
  
  gsn_panel(wks, plot, (/ 2, 2 /), resP)
  
delete(wks)
