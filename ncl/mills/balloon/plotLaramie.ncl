load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CODE_PATH/ncl/lib/boxplot_vertical.ncl"

begin
  datapath=getenv("DATA_PATH")
  imagepath=getenv("IMG_PATH")
  startyear=1999
  endyear=2008
  infile=datapath+"/Deshler2014/US_Laramie_41N_105W_"+startyear+"-"+endyear+".nc"

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


  pdffile=imagepath+"/US_Laramie_41N_105W_"+startyear+"-"+endyear
  print("rm "+pdffile+".pdf")
  rm=systemfunc("rm "+pdffile+".pdf")
  wks = gsn_open_wks("pdf",pdffile)
;**********************************************
; resources for plot background
;**********************************************
  res            = True                         ; plot mods desired
;  res@tmYLLabels = (/"Control","-2Xna","2Xna"/) ; labels for each box
  res@tiMainString = "Laramie (41N, 105W), "+startyear+"-"+endyear
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
  boxColors = (/"black","red","green4","blue","orchid","cyan4"/)  	; Color of box(es)
  opti@boxColors = boxColors(0:nbins-1)
;***********************************************
  plot = boxplot_vertical(wks,xval,alt,opti,res,llres)	; All 3 options used...
  
  print("drawing plot...")
  draw(plot)                                     ; box plot does not call
  print("framing wks...")
  frame(wks)                                    ; these for you

  print("deleting wks...")
  delete(wks)
  print(pdffile+".pdf")
end 
