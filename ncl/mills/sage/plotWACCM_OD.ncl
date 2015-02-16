load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

caseid=getenv("caseid")
if (ismissing(caseid)) then
  caseid="f.e13.FSDW5.f19_f19.1990-2011.ASv1.P17Tg.16-18km.wCH.01"
  print("missing environment variable caseid, using default")
end if
print("caseid="+caseid)

rootpath="/glade/scratch/mmills/archive/"
subdir="/atm/proc/h0zm/"

infiles=systemfunc("ls "+rootpath+caseid+subdir+"*.nc")
history=addfiles(infiles,"r")
ListSetType(history, "cat")

date=addfiles_GetVar(history,infiles,"date")
print(date)

time=yyyymmdd_to_yyyyfrac(date,0)

AODNIRst=addfiles_GetVar(history,infiles,"AODNIRst")
AODNIRst&time=time
printVarSummary(AODNIRst)
AODVISst=addfiles_GetVar(history,infiles,"AODVISst")
AODVISst&time=time
printVarSummary(AODVISst)
AODUVst=addfiles_GetVar(history,infiles,"AODUVst")
AODUVst&time=time
printVarSummary(AODUVst)

a = systemfunc("mkdir -p $IMAGE_PATH/sage")
pdffile = "$IMAGE_PATH/sage/WACCM_OD."+caseid
wks = gsn_open_wks("pdf", pdffile)

;gsn_define_colormap(wks,"rainbow")

res                     = True          ; plot mods desired
;res@trXMinF=tmin
res@trXMinF=1991.0
res@trXMaxF=1995.0
res@gsnPaperOrientation = "landscape"
res@vpWidthF            = 1.2           ; change aspect ratio of plot
res@vpHeightF           = 0.55
res@cnFillOn            = True          ; turn on color
res@gsnSpreadColors     = True          ; use full range of colormap
res@gsnMaximize         = True          ; expand plot size to maximum
res@cnLevelSelectionMode = "ExplicitLevels" 
res@cnLevels = 10.^fspan(-.5,-3.,11)

plot = gsn_csm_lat_time(wks,AODNIRst(lat|:, time|:),res)
delete(plot)

plot = gsn_csm_lat_time(wks,AODVISst(lat|:, time|:),res)
delete(plot)

plot = gsn_csm_lat_time(wks,AODUVst(lat|:, time|:),res)
delete(plot)

