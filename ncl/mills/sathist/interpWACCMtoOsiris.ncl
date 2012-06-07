load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

waccmRun="f40.Nabro.2deg.wcm.sd.carma.sulf.006"
rootdir="/Volumes/Data/Models/cesm/run/"+waccmRun"/"

hs=addfile(rootdir+"hs/"+waccmRun,"r")

osiris=addfile("$DATA_PATH/Osiris/osiris_hs_nabro_20110601_20110822.nc","r")

lev_m=osiris->lev ; Osiris levels (lev=45) in meters. We will interpolate WACCM output to these levels

z3_waccm=hs->Z3 ; WACCM geopotential height (ncol=16402, lev=88) in meters.

ext_waccm=hs->PURSULEX750 ; WACCM extinction (ncol=16402, lev=88) in /km. 

; Flip WACCM levels and interpolate to Osiris levels.
ext_waccm_osiris_levs = linint1_Wrap (z3_waccm(:,::-1), ext_waccm(:,::-1), False, lev_m, 0)

outfile=addfile(rootdir+"hs/"+waccmRun+".cam.hs.OsirisLevs.nc","c")
outfile->PURSULEX750=ext_waccm_osiris_levs

