load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CODE_PATH/ncl/lib/boxplot_vertical.ncl"

begin
  datapath=getenv("DATA_PATH")
  imagepath=getenv("IMG_PATH")
  rangestr=getenv("rangestr")
  caseid=getenv("caseid")
  
  if (ismissing(rangestr)) then
    rangestr="1999-2008"
  end if  
    
  print("Laramie range: "+rangestr)
  
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
    modelrange=getenv("modelrange")
    if ismissing(modelrange) then
      modelrange=rangestr
      pdffile=pdffile+"."+caseid
    else
      pdffile=pdffile+"."+caseid+"_"+modelrange
    end if
    print("model using year(s) "+rangestr+" from case "+caseid)
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
;  res@tiMainOffsetYF    = 0.02
;  res@tiMainFontHeightF = 0.035
  res@tiYAxisString = "Altitude (km)"
  res@tiXAxisString = "Number concentration (cm-3)"
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
  
  if ismissing(caseid) then
    plot = boxplot_vertical(wks,xval,alt,opti,res,llres)	; All 3 options used...
    print("drawing plot...")
    draw(plot)                                     ; box plot does not call
    print("framing wks...")
    frame(wks)                                    ; these for you    
  else
    plot = new(4, graphic)    
    plot(0) = 	boxplot_vertical(wks,xval,alt,opti,res,llres)
    res@tiMainString = "Model mode 1, accumulation"
    plot(1) = 	boxplot_vertical(wks,xval,alt,opti,res,llres)
    res@tiMainString = "Model mode 2, Aitken"
    plot(2) = 	boxplot_vertical(wks,xval,alt,opti,res,llres)
    res@tiMainString = "Model mode 3, coarse"
    plot(3) = 	boxplot_vertical(wks,xval,alt,opti,res,llres)
    modelfile=caseid+".US_Laramie_41N_105W_"+modelrange+".nc"
    rootpath = "/glade/scratch/mmills/archive/"

    exists =isfilepresent(rootpath+caseid+"/atm/proc/"+modelfile)
    if (.not.exists) then
      print("Missing: "+rootpath+caseid+"/atm/proc/"+modelfile)
      rootpath="/glade/p/cesm/wawg_dev/mmills/archive/"
    end if

    exists =isfilepresent(rootpath+caseid+"/atm/proc/"+modelfile)
    if (.not.exists) then
      print("Missing: "+rootpath+caseid+"/atm/proc/"+modelfile)
      print("drawing plot...")
      draw(plot)                                     ; box plot does not call
      print("framing wks...")
      frame(wks)                                    ; these for you  
    else
      print("opening "+rootpath+caseid+"/atm/proc/"+modelfile)
      ncid=addfile(rootpath+caseid+"/atm/proc/"+modelfile,"r")
      ndbin=ncid->ndbin
      nd1bin=ncid->nd1bin
      nd2bin=ncid->nd2bin
      nd3bin=ncid->nd3bin
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
      overlay(plot(0),gsn_csm_xy(wks,ndbin(bin|:,lev|:),modalt,res))
      overlay(plot(1),gsn_csm_xy(wks,nd1bin(bin|:,lev|:),modalt,res))
      overlay(plot(2),gsn_csm_xy(wks,nd2bin(bin|:,lev|:),modalt,res))
      overlay(plot(3),gsn_csm_xy(wks,nd3bin(bin|:,lev|:),modalt,res))
;
; Create a string to put at the top, using a slightly
; smaller font height than the axis titles.
;
;      txres = True
;      txres@txFontColor="black"
;      txres@txFontHeightF = 0.015  ; font smaller. default big
;      txidc = gsn_create_text(wks, caseid+", "+modelrange, txres)
;      amres                  = True
;      amres@amOrthogonalPosF = -0.53   ; Top of plot plus a little extra                                     ; to stay out of the tickmarks.
;      amres@amJust           = "BottomCenter"
;      amres@amParallelPosF   = 0.0   ; Centered
;      annoidl = gsn_add_annotation(plot(0), txidc, amres)
;      txres@txFontHeightF = 0.02  ; font smaller. default big
;      txidc = gsn_create_text(wks, "model mode 1, accumulation", txres)
;      annoidl1 = gsn_add_annotation(plot(1), txidc, amres)
;      txidc = gsn_create_text(wks, "model mode 2, Aitken", txres)
;      annoidl2 = gsn_add_annotation(plot(2), txidc, amres)
;      txidc = gsn_create_text(wks, "model mode 3, coarse", txres)
;      annoidl3 = gsn_add_annotation(plot(3), txidc, amres)
;      print((/ndbin&bin/))
;      print(sprintf("%3.1f",modalt)+" "+sprintf("%3.2e",ndbin(:,0))+" "+sprintf("%3.2e",ndbin(:,1))+" "+\
;        sprintf("%3.2e",ndbin(:,2))+" "+sprintf("%3.2e",ndbin(:,3))+" "+sprintf("%3.2e",ndbin(:,4))+" "+\
;        sprintf("%3.2e",ndbin(:,5))+" "+sprintf("%3.2e",ndbin(:,6))+" "+sprintf("%3.2e",ndbin(:,7)))

      resP                       = True
      resP@gsnMaximize           = True            ; expand plot size to maximum
      resP@gsnPaperOrientation   = "portrait"
      resP@txString   = caseid+", "+modelrange

      gsn_panel(wks, plot, (/ 2, 2 /), resP)
      
    end if
  end if

  print("deleting wks...")
  delete(wks)
  print(pdffile+".pdf")
end 
