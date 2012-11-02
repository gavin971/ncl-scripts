load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("plotVars")
procedure plotVars(path:string,cases:string,varname:string,times:integer,lat:float,lon:float)

begin
  ncases=dimsizes(cases)
  data=new((/ncases,times/),"float")
  do i = 0,ncases-1
    print(i+" "+cases(i))
    files=systemfunc("ls "+path+"/"+cases(i)+"/h0/*.nc")
    f=addfiles(files(0:times-1),"r")
    var=f[:]->$varname$
    data(i,:)=var(:,{lat},{lon})
    time=f[:]->time
    delete(files)
    delete(f)
    delete(var)
  end do
  
  pdffile=cases(0)+"."+varname
  wks = gsn_open_wks ("pdf","$IMG_PATH/"+pdffile)
  res=True
  res@trXMaxF=time(times-1)
  res@xyMonoLineColor   = False
  res@xyLineColors   = (/"Black","Red","Blue"/)
  res@xyMonoDashPattern   = True
  plot=gsn_csm_xy(wks,time,data,res)
  delete(wks)
  
end

path="/glade/proj3/CESM/wawg/runs"
cases=(/"b40_1850_2d_r07c5cn_160jp","b.e11.B1850W5CN.f19_g16.007a","b.e11.B1850W5CN.f19_g16.007b"/)
varname="PBLH"
times=7
lat=25.5
lon=32.5

plotVars(path,cases,varname,times,lat,lon)
