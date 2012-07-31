load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  

filename="$DATA_PATH/HADCRU/HadCRUT3v.nc"
hadcru=addfile(filename,"r")

temp=hadcru->temp(:,0,:,:)

d1=0.0
d2=365.0*20.
temp_PI  =dim_avg_Wrap(temp(latitude|:,longitude|:,{t|d1:d2}))

d1=(1985-1850)*365.0
d2=(2005-1850)*365.0
temp_20th=dim_avg_Wrap(temp(latitude|:,longitude|:,{t|d1:d2}))

temp_diff = temp_20th 
temp_diff = (/temp_20th - temp_PI/)


  wks = gsn_open_wks("eps","$IMG_PATH/hadcrut3")                  ; open a ps file
  plot = gsn_csm_contour_map_ce(wks,temp_diff,False)     ; create a default plot
