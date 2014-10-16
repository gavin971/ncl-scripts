filename="/glade/scratch/mmills/archive/b.e13.B20TRW5TCN.f19_g16.igwsProgAer.001/atm/proc/h0gm/b.e13.B20TRW5TCN.f19_g16.igwsProgAer.001.cam.h0.globavg.nc"

a=addfile(filename,"r")

CH4_CHML=a->CH4_CHML

r_CH4_OH=a->r_CH4_OH
r_CL_CH4=a->r_CL_CH4
r_F_CH4=a->r_F_CH4
r_O1D_CH4a=a->r_O1D_CH4a
r_O1D_CH4b=a->r_O1D_CH4b
r_O1D_CH4c=a->r_O1D_CH4c
r_jch4_a=a->r_jch4_a
r_jch4_b=a->r_jch4_b

z=a->Z3

z=z/1000.

ds=dimsizes(CH4_CHML)
print(ds)

t=0 

; Open a workstation
image_path=getenv("IMAGE_PATH")
pdffile=image_path+"/CH4lossVertProf"
print("Creating PDF file "+pdffile+".pdf")
wks   = gsn_open_wks ("pdf", pdffile)

res                       = True
res@gsnMaximize = True
res@vpWidthF              = 0.6
res@vpHeightF             = 0.8
res@gsnDraw               = False
res@gsnFrame              = False
res@trXLog                = True
res@trXMinF               = 1.e-3
res@trXMaxF               = 2e5
res@trYMinF               = 0.
res@trYMaxF               = 110.
res@xyLineThicknesses = (/ 3 /)
res@tiXAxisString = CH4_CHML@long_name+" ("+CH4_CHML@units+")"
res@tiYAxisString = "Altitude (km)"

res@xyLineColors          = (/ "black" /)  
plot=gsn_csm_xy(wks, CH4_CHML(t,:), z(t,:), res)

res@xyLineColors          = (/ "red" /)  
overlay(plot,gsn_csm_xy(wks, r_CH4_OH(t,:), z(t,:), res))

res@xyLineColors          = (/ "orange" /)  
overlay(plot,gsn_csm_xy(wks, r_CL_CH4(t,:), z(t,:), res))

res@xyLineColors          = (/ "yellow" /)  
overlay(plot,gsn_csm_xy(wks, r_F_CH4(t,:), z(t,:), res))

res@xyLineColors          = (/ "greenyellow" /)  
overlay(plot,gsn_csm_xy(wks, r_O1D_CH4a(t,:), z(t,:), res))

res@xyLineColors          = (/ "green" /)  
overlay(plot,gsn_csm_xy(wks, r_O1D_CH4b(t,:), z(t,:), res))

res@xyLineColors          = (/ "darkgreen" /)  
overlay(plot,gsn_csm_xy(wks, r_O1D_CH4c(t,:), z(t,:), res))

res@xyLineColors          = (/ "blue" /)  
overlay(plot,gsn_csm_xy(wks, r_jch4_a(t,:), z(t,:), res))

res@xyLineColors          = (/ "purple" /)  
overlay(plot,gsn_csm_xy(wks, r_jch4_b(t,:), z(t,:), res))

res@xyLineThicknesses = (/ 5 /)
res@xyLineColors          = (/ "black" /)  
overlay(plot,gsn_csm_xy(wks, CH4_CHML(t,:), z(t,:), res))

  xpos=2.0e-3
  ypos=50.0
  ydif=5.0
  txres = True
  txres@txFontColor="black"
  txres@txFontHeightF = 0.02                     ; font smaller. default big
  txres@txJust="CenterLeft"
  dum=gsn_add_text(wks,plot(0),"CH4 loss",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="red"
  dum=gsn_add_text(wks,plot(0),"CH4 + OH -> CH3O2 + H2O",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="orange"
  dum=gsn_add_text(wks,plot(0),"Cl + CH4 -> CH3O2 + HCl",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="yellow"
  dum=gsn_add_text(wks,plot(0),"F + CH4 -> HF + CH3O2",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="greenyellow"
  dum=gsn_add_text(wks,plot(0),"O1D + CH4 -> CH3O2 + OH",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="green"
  dum=gsn_add_text(wks,plot(0),"O1D + CH4 -> CH2O + H + HO2",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="darkgreen"
  dum=gsn_add_text(wks,plot(0),"O1D + CH4 -> CH2O + H2",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="blue"
  dum=gsn_add_text(wks,plot(0),"CH4 + hv -> H + CH3O2",xpos,ypos,txres)
  
  ypos=ypos-ydif
  txres@txFontColor="purple"
  dum=gsn_add_text(wks,plot(0),"CH4 + hv -> 1.44*H2 + ...",xpos,ypos,txres)

draw(plot)
frame(wks)

delete(wks)

print("evince "+pdffile+".pdf")
evincer=systemfunc("evince " + pdffile + ".pdf &")
print(evincer)

