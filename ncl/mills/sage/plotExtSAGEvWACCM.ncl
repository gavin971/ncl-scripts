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
  
ostype=getenv("OSTYPE")
host=getenv("HOSTNAME")

rootpath="/glade/scratch/mmills/archive/"
subdir="/atm/proc/h0zm/"

exists =isfilepresent(rootpath+caseid+subdir)
if (.not.exists) then
  print("Missing: "+rootpath+caseid+subdir)
  rootpath="/glade/p/cesm/wawg_dev/mmills/archive/"
end if

exists =isfilepresent(rootpath+caseid+subdir)
if (.not.exists) then
  print("Missing: "+rootpath+caseid+subdir)
  print("FATAL: check caseid="+caseid)
  exit
end if
print("Rootpath: "+rootpath+caseid+subdir)

sagepath="/glade/p/work/mmills/data/SAGE_II_v7.00/CCMI_Filled_Data_V2/netCDF/"
if (ostype.eq."darwin") then
  sagepath="/Volumes/Data/Data/SAGE2/SAGE_II_v7.00/CCMI_Filled_Data_V2/netCDF/"
end if


b=addfile(sagepath+"CCMI_Filled_Data_V2.nc","r")
ext1020s=b->ext1020
ext525s=b->ext525
ext452s=b->ext452
ext386s=b->ext386

sageLat=ext1020s&lat
nLat=dimsizes(sageLat)
sageAlt=ext1020s&alt
nAlt=dimsizes(sageAlt)

sagedate=ext1020s&time
sage1991=ind(floattoint(sagedate).eq.1991)
sageD1=sage1991(4)
print(sageD1+" "+sagedate(sageD1))

infiles=systemfunc("ls "+rootpath+caseid+subdir+"*.nc")
history=addfiles(infiles,"r")
ListSetType(history, "cat")

waccmdate=addfiles_GetVar(history,infiles,"date")
waccmD1=ind(waccmdate.eq.19910601)
print(waccmD1+" "+waccmdate(waccmD1))

ext550w=addfiles_GetVar(history,infiles,"EXTINCT")
ext1020w=addfiles_GetVar(history,infiles,"EXTINCTNIR")
waccmLat=ext550w&lat
ext550w=ext550w*1e3
ext550w@units="km-1"
ext1020w=ext1020w*1e3
ext1020w@units="km-1"
;print(ext1020w&lev@long_name)
ext1020w&lev@long_name="Pressure"
ext550w&lev@long_name="Pressure"

;nYears=4
;nMonths=nYears*12
nMonths=8
date=sagedate(sageD1:sageD1+nMonths-1)

a = systemfunc("mkdir -p $IMAGE_PATH/sage")
pdffile = "$IMAGE_PATH/sage/ExtSAGEvWACCM."+caseid
wks = gsn_open_wks("pdf", pdffile)

resP=True
resP@gsnPaperOrientation = "landscape"
resP@gsnMaximize         = True          ; expand plot size to maximum
resP@gsnPanelRight =0.95
;resP@gsnPanelFigureStrings =(/"SAGE 1020nm","WACCM 1020nm","SAGE 525nm","WACCM 550nm"/)

res=True
res@cnFillOn            = True          ; turn on color
res@gsnSpreadColors     = True          ; use full range of colormap
res@gsnDraw = False
res@gsnFrame = False
res@vpWidthF            = 1.1           ; change aspect ratio of plot
res@vpHeightF           = 0.7
;res@gsnMaximize         = True          ; expand plot size to maximum
res@cnLevelSelectionMode = "ExplicitLevels" 
res@cnLevels = 10.^fspan(-.75,-3.5,12)
res@gsnLeftString=" "
res@gsnRightString=" "
res@gsnCenterStringFontHeightF=0.04

resW=res
resW@trYMinF=10
resW@trYMaxF=270
resW@trXMinF=-77
resW@trXMaxF=77
;resW@tiYAxisString="Pressure (hPa)"

res@trYMaxF=31
res@trYMinF=10

months=(/"January", "February", "March", "April", "May", "June", "July", \
       "August", "September", "October", "November", "December" /)
plot = new(4, graphic)
do year=1991,1991
  do month=5,12
    d=(year-1991)*12+month-5
    wd=waccmD1+d
    sd=sageD1+d
    print(d+" "+wd+" "+sd)
    
    res@gsnCenterString="SAGE "+ext1020s@long_name+" (/km)"
    plot(0)  = gsn_csm_contour(wks, ext1020s(sd,:,:), res ) 
    res@gsnCenterString="SAGE "+ext525s@long_name+" (/km)"
    plot(2)  = gsn_csm_contour(wks, ext525s(sd,:,:), res ) 
    resW@gsnCenterString="WACCM "+ext1020w@long_name+" (/km)"
    plot(1)  = gsn_csm_pres_hgt(wks, ext1020w(wd,:,:), resW ) 
    resW@gsnCenterString="WACCM "+ext550w@long_name+" (/km)"
    plot(3)  = gsn_csm_pres_hgt(wks, ext550w(wd,:,:), resW ) 
;    txres = True
;    txres@txFontColor="black"
;    txres@txFontHeightF = 0.05                     ; font smaller. default big
;    xpos=0
;    ypos=30.0
    datestring=months(month-1)+" "+year
   ; dum = gsn_add_text(wks,plot(0),datestring,xpos,ypos,txres) 
    resP@gsnPanelFigureStrings =(/"",datestring,"",""/)
    resP@amJust="TopLeft"
    resP@gsnPanelFigureStringsFontHeightF = 0.018
    
    gsn_panel(wks, plot, (/ 2, 2 /), resP)
;    frame(wks)
  end do
end do

delete(wks)
