load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  

rootPath="/data/ccsm4/"
filename1=rootPath+"b40.1850.track1.1deg.006/climo/1271-1300/b40.1850.track1.1deg.006_ANN_climo.nc"
filename2=rootPath+"b40.20th.track1.1deg.005/climo/1986-2004/b40.20th.track1.1deg.005_ANN_climo.nc"

filePI=addfile(filename1,"r")
file20=addfile(filename2,"r")

tempPI=filePI->TS
temp20=file20->TS


temp_diff = temp20(0,:,:)
temp_diff = (/temp20(0,:,:) - tempPI(0,:,:)/)


wks = gsn_open_wks("eps","$IMG_PATH/ccsm4TS")                  ; open a ps file
plot = gsn_csm_contour_map_ce(wks,temp_diff,False)     ; create a default plot
