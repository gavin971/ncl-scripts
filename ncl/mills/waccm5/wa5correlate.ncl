load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

lat=25.5
lon=32.5

rootdir = "/glade/proj3/CESM/wawg/runs/"
waccm5_case="b.e11.B1850W5CN.f19_g16.007"
cam5_case="b40_1850_2d_r07c5cn_160jp"

print("reading WACCM files")
infiles_wcm = systemfunc("ls "+rootdir+waccm5_case+"/h0/"+waccm5_case+".cam.h0.*-07.nc")
history_wcm = addfiles(infiles_wcm, "r")
ListSetType(history_wcm, "cat")
U10_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "U10")
DSTSFMBL_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "DSTSFMBL")
TS_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "TS")
U_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "U")
V_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "V")
date_wcm = addfiles_GetVar(history_wcm, infiles_wcm, "date")

print("reading CAM files")
infiles_cam = systemfunc("ls "+rootdir+cam5_case+"/h0/"+cam5_case+".cam2.h0.*-07.nc")
history_cam = addfiles(infiles_cam, "r")
ListSetType(history_cam, "cat")
U10_cam = addfiles_GetVar(history_cam, infiles_cam, "U10")
DSTSFMBL_cam = addfiles_GetVar(history_cam, infiles_cam, "DSTSFMBL")
TS_cam = addfiles_GetVar(history_cam, infiles_cam, "TS")
U_cam = addfiles_GetVar(history_cam, infiles_cam, "U")
V_cam = addfiles_GetVar(history_cam, infiles_cam, "V")

  ; Plot the correlation of the data.
  pdfpath="$IMG_PATH/"+waccm5_case+"/"
  a = systemfunc("mkdir -p "+pdfpath)
  pdffile=pdfpath+waccm5_case+".corr."+cam5_case+".lat"+sprintf("%04.1f",lat)+".lon"+sprintf("%04.1f",lon)
  print("Creating file "+pdffile+".pdf")
  wks  = gsn_open_wks("pdf", pdffile)    ; open a ncgm file
  gsn_define_colormap(wks, "rainbow+white+gray")  ; choose a colormap

  res                   = True
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@tiMainString      = "Lat: "+sprintf("%4.1f",lat)+" Lon: "+sprintf("%4.1f",lon)
  res@xyMarkLineModes   = "Markers"                ; choose which have markers
  res@xyMarkers         =  16                      ; choose type of marker  
  res@xyMarkerColor     = "red"                    ; Marker color
  res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
;  res@tiXAxisString     = "WACCM/CARMA model"
;  res@tiYAxisString     = "Osiris observations"
  res@trXMinF           = 0.0
  res@trYMinF           = 0.0
  res@trXMaxF           = 5.0
  res@trYMaxF           = 1.5e-6
  
  print("Plotting scatterplot...")
  
  plot = gsn_csm_xy(wks, U10_wcm(:,{lat},{lon}), DSTSFMBL_wcm(:,{lat},{lon}), res)       ; scatterplot
  res@xyMarkerColor     = "blue"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, U10_cam(:,{lat},{lon}), DSTSFMBL_cam(:,{lat},{lon}), res))       ; scatterplot
  draw(plot)
  frame(wks)
  
  res@trXMinF           = 300.0
  res@trXMaxF           = 318.0
  res@xyMarkerColor     = "red"                    ; Marker color
  plot = gsn_csm_xy(wks, TS_wcm(:,{lat},{lon}), DSTSFMBL_wcm(:,{lat},{lon}), res)       ; scatterplot
  res@xyMarkerColor     = "blue"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, TS_cam(:,{lat},{lon}), DSTSFMBL_cam(:,{lat},{lon}), res))       ; scatterplot
  draw(plot)
  frame(wks)
    
  res@trYMinF           = 2.5
  res@trYMaxF           = 5.0
  res@xyMarkerColor     = "red"                    ; Marker color
  plot = gsn_csm_xy(wks, TS_wcm(:,{lat},{lon}), U10_wcm(:,{lat},{lon}), res)       ; scatterplot
  res@xyMarkerColor     = "blue"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, TS_cam(:,{lat},{lon}), U10_cam(:,{lat},{lon}), res))       ; scatterplot
  draw(plot)
  frame(wks)
  
  res@trXMinF           = 0.0
  res@trXMaxF           = 3.0
  res@trYMinF           = -7.0
  res@trYMaxF           = 0.0
  res@xyMarkerColor     = "red"                    ; Marker color
  plot = gsn_csm_xy(wks, U_wcm(:,69,{lat},{lon}), V_wcm(:,69,{lat},{lon}), res)       ; scatterplot
  res@xyMarkerColor     = "blue"                    ; Marker color
  overlay(plot, gsn_csm_xy(wks, U_cam(:,29,{lat},{lon}), V_cam(:,29,{lat},{lon}), res))       ; scatterplot
  draw(plot)
  frame(wks)
  
  delete(res@trXMinF)
  delete(res@trXMaxF)
  delete(res@trYMinF)
  delete(res@trYMaxF)
  res@xyMarkerColor     = "red"                    ; Marker color
  plot = gsn_csm_xy(wks, date_wcm/10000, U10_wcm(:,{lat},{lon}), res)       ; scatterplot
  draw(plot)
  frame(wks)
  
  
  delete(wks)
  delete(plot)
  delete(res)
