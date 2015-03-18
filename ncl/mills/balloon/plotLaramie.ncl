load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CODE_PATH/ncl/lib/boxplot_vertical.ncl"

begin
  datapath=getenv("DATA_PATH")
  imagepath=getenv("IMG_PATH")
  startyear=getenv("startyear")
  endyear=getenv("endyear")
  caseid=getenv("caseid")
  
  if (ismissing(startyear)) then
    startyear=1999
  end if  
  
  if (ismissing(endyear)) then
    endyear=2008
  end if
  
  if (startyear.gt.endyear) then
    print("FATAL: startyear = "+startyear+" > endyear = "+endyear)
    exit
  end if
  
  print("startyear = "+startyear+", endyear = "+endyear)
  
  if (startyear.eq.endyear) then
    rangestr=startyear
  else
    rangestr=startyear+"-"+endyear
  end if
  
  infile=datapath+"/Deshler2014/US_Laramie_41N_105W_"+rangestr+".nc"

  print("opening "+infile)
  ncid=addfile(infile,"r")
  print("reading Nmedian...")
  Nmedian=ncid->Nmedian
  print("reading Nmin...")
  Nmin=ncid->Nmin
  print("reading Nmax...")
  Nmax=ncid->Nmax
  print("reading N25pct...")
  N25pct=ncid->N25pct
  print("reading N75pct...")
  N75pct=ncid->N75pct
  print("reading Nk...")
  Nk=ncid->Nk
  delete(ncid)
  
  diminq=dimsizes(Nmedian)
  nbins=diminq(0)
;  nbins=4

  minalt=10
  maxalt=32
  mink=ind(Nmedian&alt.eq.minalt)
  maxk=ind(Nmedian&alt.eq.maxalt)
  alt=Nmedian&alt(mink:maxk)

  nalts=dimsizes(alt)
  xval=new((/nbins,nalts,5/),"float")
  xval!0="bin"
  xval&bin=Nmedian&bin(0:nbins-1)
  xval!1="alt"
  xval&alt=alt
  minv=1e-5
  xval(0:nbins-1,:,0)=where(Nmin(0:nbins-1,mink:maxk).ge.minv,Nmin(0:nbins-1,mink:maxk),minv)
  xval(0:nbins-1,:,1)=where(N25pct(0:nbins-1,mink:maxk).ge.minv,N25pct(0:nbins-1,mink:maxk),minv)
  xval(0:nbins-1,:,2)=where(Nmedian(0:nbins-1,mink:maxk).ge.minv,Nmedian(0:nbins-1,mink:maxk),minv)
  xval(0:nbins-1,:,3)=where(N75pct(0:nbins-1,mink:maxk).ge.minv,N75pct(0:nbins-1,mink:maxk),minv)
  xval(0:nbins-1,:,4)=where(Nmax(0:nbins-1,mink:maxk).ge.minv,Nmax(0:nbins-1,mink:maxk),minv)

  print("mkdir "+imagepath+"/balloon")
  cmd=systemfunc("mkdir "+imagepath+"/balloon")
  pdffile=imagepath+"/balloon/US_Laramie_41N_105W_"+rangestr
  if .not.ismissing(caseid) then
    pdffile=pdffile+"."+caseid
  end if
  print("rm "+pdffile+".pdf")
  rm=systemfunc("rm "+pdffile+".pdf")
  wks = gsn_open_wks("pdf",pdffile)
;**********************************************
; resources for plot background
;**********************************************
  res            = True                         ; plot mods desired
;  res@tmYLLabels = (/"Control","-2Xna","2Xna"/) ; labels for each box
  res@tiMainString = "Laramie (41N, 105W), "+rangestr
;  res@trXLog=True
;**********************************************
; resources for polylines that draws the boxes
;**********************************************  
  llres                   = True			
  llres@gsLineThicknessF  = 2.5                 ; line thickness 
;**********************************************
; resources that control color and width of boxes
;**********************************************  
  opti          = True			
  opti@boxWidth = .25				; Width of box (y units)
;  boxColors = (/"black","red","green4","blue","orchid","orange","cyan4","grey"/)  	; Color of box(es)
  boxColors = (/"black","red","green4","blue","orchid","transparent","cyan4","transparent"/)  	; Color of box(es)
  opti@boxColors = boxColors(0:nbins-1)
;***********************************************
  plot = boxplot_vertical(wks,xval,alt,opti,res,llres)	; All 3 options used...
  
  
  if .not.ismissing(caseid) then
    modelfile=caseid+".US_Laramie_41N_105W_"+rangestr+".nc"
    rootpath = "/glade/scratch/mmills/archive/"

    exists =isfilepresent(rootpath+caseid+"/atm/proc/"+modelfile)
    if (.not.exists) then
      print("Missing: "+rootpath+caseid+"/atm/proc/"+modelfile)
      rootpath="/glade/p/cesm/wawg_dev/mmills/archive/"
    end if

    exists =isfilepresent(rootpath+caseid+"/atm/proc/"+modelfile)
    if (.not.exists) then
      print("Missing: "+rootpath+caseid+"/atm/proc/"+modelfile)
    else
      print("opening "+rootpath+caseid+"/atm/proc/"+modelfile)
      ncid=addfile(rootpath+caseid+"/atm/proc/"+modelfile,"r")
      ndbin=ncid->ndbin
      Z3=ncid->Z3
      modalt=Z3/1e3
      print(xval&bin)
      print(ndbin&bin)
      print("overlaying plot...")
      res@xyMonoDashPattern=True
      res@xyLineThicknessF=3.0
      res@xyLineColors=boxColors
      res@xyMonoLineColor=False
      res@gsnDraw=False
      res@gsnFrame=False
      overlay(plot,gsn_csm_xy(wks,ndbin(bin|:,lev|:),modalt,res))
    end if
  end if

  print("drawing plot...")
  draw(plot)                                     ; box plot does not call
  print("framing wks...")
  frame(wks)                                    ; these for you

  print("deleting wks...")
  delete(wks)
  print(pdffile+".pdf")
end 
