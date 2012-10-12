load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

lat=25.5
lon=32.5

rootdir = "/glade/proj3/CESM/wawg/runs/"
waccm5_case="b.e11.B1850W5CN.f19_g16.009"

print("ls "+rootdir+waccm5_case+"/h0/*.h0.*.nc")
infiles_wcm = systemfunc("ls "+rootdir+waccm5_case+"/h0/*.h0.*.nc")
print("addfiles(infiles_wcm)")
history_wcm = addfiles(infiles_wcm, "r")
print("ListSetType(history_wcm)")
ListSetType(history_wcm, "cat")
print("reading U10")
U10_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "U10")
print("reading DSTSFMBL")
DSTSFMBL_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "DSTSFMBL")
print("reading TS")
TS_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "TS")
print("reading U")
U_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "U")
print("reading V")
V_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "V")
;date_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "date")
print("reading time")
time_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "time")

  ; Plot the correlation of the data.
  pdfpath="$IMG_PATH/"+waccm5_case+"/"
  a = systemfunc("mkdir -p "+pdfpath)
  pdffile=pdfpath+waccm5_case+".inst.lat"+sprintf("%04.1f",lat)+".lon"+sprintf("%04.1f",lon)
  print("Creating file "+pdffile+".pdf")
  wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow+white+gray")  ; choose a colormap

  res                   = True
;  res@gsnDraw           = False
;  res@gsnFrame          = False
  res@tiMainString      = "Lat: "+sprintf("%4.1f",lat)+" Lon: "+sprintf("%4.1f",lon)
  res@xyMarkLineModes   = "Markers"                ; choose which have markers
  res@xyMarkers         =  16                      ; choose type of marker  
  res@xyMarkerColor     = "red"                    ; Marker color
  res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
;  res@tiXAxisString     = "WACCM/CARMA model"
;  res@tiYAxisString     = "Osiris observations"
;  res@trXMinF           = 0.0
;  res@trYMinF           = 0.0
;  res@trXMaxF           = 5.0
;  res@trYMaxF           = 1.5e-6
  
  print("Plotting...")
  plot = gsn_csm_xy(wks, time_wcm, U10_wcm(:,{lat},{lon}), res)       ; scatterplot
  plot = gsn_csm_xy(wks, time_wcm, U_wcm(:,69,{lat},{lon}), res)       ; scatterplot
  plot = gsn_csm_xy(wks, time_wcm, V_wcm(:,69,{lat},{lon}), res)       ; scatterplot
  plot = gsn_csm_xy(wks, time_wcm, TS_wcm(:,{lat},{lon}), res)       ; scatterplot
  plot = gsn_csm_xy(wks, time_wcm, DSTSFMBL_wcm(:,{lat},{lon}), res)       ; scatterplot
  
  delete(wks)
  delete(plot)
  delete(res)
