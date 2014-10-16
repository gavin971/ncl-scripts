load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

run_name = getenv("caseid")

if (ismissing(run_name)) then
  print("set environment variable caseid before running")
else

  year0 = 1
;  nYears=10-year0+1

  ostype=getenv("OSTYPE")
  host=getenv("HOSTNAME")
  
  subdir=getenv("subdir")
  if (ismissing(subdir)) then
    subdir="atm/hist"
  end if

  rootpath=getenv("rootpath")
  if (ismissing(rootpath)) then
    if (ostype.eq."darwin") then
      rootpath="/Volumes/MillsPassport/cesm/"
    else
      if(isStrSubset(host,"cmip5"))
        rootpath="/data/waccm/"
      else
	if(isStrSubset(host,"lfe")) ; lou (NAS HECC)
          rootpath="~/csm/"
	else
          if(isStrSubset(host,"pfe"))
            rootpath="/nobackupp8/mjmills2/archive/"
          else
            rootpath="/glade/scratch/mmills/archive/"
            exists =isfilepresent(rootpath+run_name+"/"+subdir+"/")
            if (.not.exists) then
              rootpath="/glade/p/cesm/wawg/runs/cesm1_2/"
              delete(exists)
              exists =isfilepresent(rootpath+run_name+"/"+subdir+"/")
              if (.not.exists) then
        	rootpath="/glade/p/cesm/wawg/runs/ccsm4_0/"
              end if
            end if  
          end if  
        end if  
      end if
    end if
  end if

  in_dir =rootpath+run_name+"/"+subdir+"/"
  print("in_dir="+in_dir)
  save_dir =rootpath+run_name+"/atm/proc/"
  print("mkdir -p "+save_dir)
  a = systemfunc("mkdir -p "+save_dir)
  savefile=save_dir+run_name+"_RESTOM.nc"
  print("savefile="+savefile)

  image_path=getenv("IMAGE_PATH")
  pdffile=image_path+"/"+run_name+"_RESTOM"
  print("Creating PDF file "+pdffile+".pdf")
  wks   = gsn_open_wks ("pdf", pdffile)

  mon_wgt = (/0.08493150770664215,0.07671232521533966,0.08493150770664215,\
              0.08219178020954132,0.08493150770664215,0.08219178020954132,\
              0.08493150770664215,0.08493150770664215,0.08219178020954132,\
              0.08493150770664215,0.08219178020954132,0.08493150770664215/)

  infiles = systemfunc("ls -1 "+in_dir+"*h0.*.nc")
  nYears=dimsizes(infiles)/12
  print("nYears="+nYears)

  nYears_a=0
  overwrite=getenv("overwrite")
  doTropFlux = False
  
  if (isfilepresent(savefile).and.ismissing(overwrite)) then

    print("reading "+savefile)
    b = addfile(savefile,"r")
    year_a=b->year
    restom_ann_w_a=b->restom_ann_w
    ts_ann_w_a=b->ts_ann_w
    ts_land_ann_w_a=b->ts_land_ann_w
    sst_ann_w_a=b->sst_ann_w

    ; query if clearsky restom exists
    if (isfilevar(b,"restom_ann_w_c")) then
      print("found restom_ann_w_c")
      restom_ann_w_c_a=b->restom_ann_w_c
      flnt_ann_w_a=b->flnt_ann_w
      flntc_ann_w_a=b->flntc_ann_w
      fsnt_ann_w_a=b->fsnt_ann_w
      fsntc_ann_w_a=b->fsntc_ann_w
      nYears_a = dimsizes(year_a)
    end if
    
    if (isfilevar(b,"flnr_ann_w")) then
      print("found flnr_ann_w")
      flnr_ann_w_a=b->flnr_ann_w
      fsnr_ann_w_a=b->fsnr_ann_w
      doTropFlux = True
    end if
    
  end if

  if ((nYears_a.ge.nYears)) then ; use only savefile if it has all history years

    year=year_a
    restom_ann_w=restom_ann_w_a
    ts_ann_w=ts_ann_w_a
    ts_land_ann_w=ts_land_ann_w_a
    sst_ann_w=sst_ann_w_a  

  else

    i=nYears_a*12
    fin=addfile(infiles(i),"r")
    print((/infiles(i)/))
  ;  print(fin->date)
    fsnt = fin->FSNT
    flnt = fin->FLNT
    fsntc = fin->FSNTC
    flntc = fin->FLNTC
    ts = fin->TS
    landfrac = fin->LANDFRAC
    year_in = fin->date ; get record coordinate
    year_in = year_in/10000
  ;  print(year_in)
    if (isfilevar(fin,"FLNR")) then
      fsnr = fin->FSNR
      flnr = fin->FLNR
      doTropFlux = True
    end if

    do i=nYears_a*12+1,dimsizes(infiles)-1,1
      print(i+" "+(/infiles(i)/))
      fin=addfile(infiles(i),"r")
  ;    print(fin->date)

      temp=fsnt
      delete(fsnt)
      fsnt = array_append_record(temp,fin->FSNT,0)
      delete(temp)

      temp=flnt
      delete(flnt)
      flnt = array_append_record(temp,fin->FLNT,0)
      delete(temp)

      temp=fsntc
      delete(fsntc)
      fsntc = array_append_record(temp,fin->FSNTC,0)
      delete(temp)

      temp=flntc
      delete(flntc)
      flntc = array_append_record(temp,fin->FLNTC,0)
      delete(temp)

      temp=ts
      delete(ts)
      ts = array_append_record(temp,fin->TS,0)
      delete(temp)

      temp=landfrac
      delete(landfrac)
      landfrac = array_append_record(temp,fin->LANDFRAC,0)
      delete(temp)

      temp=year_in
      delete(year_in)
  ;    print(i+" "+temp+" "+fin->date+" "+(fin->date)/10000)
      year_in = array_append_record(temp,(fin->date)/10000,0)
      delete(temp)
      
      if (doTropFlux) then
        temp=fsnr
        delete(fsnr)
        fsnr = array_append_record(temp,fin->FSNR,0)
        delete(temp)
        
        temp=flnr
        delete(flnr)
        flnr = array_append_record(temp,fin->FLNR,0)
        delete(temp)
      end if

    end do
    ;    fid=addfile(in_file,"r")
    ;    data = fid->$variable$
    ;inptr = addfiles(infiles,"r")
    ;fsnt    = addfiles_GetVar(inptr,infiles,"FSNT")
    ;flnt    = addfiles_GetVar(inptr,infiles,"FLNT")

    ;ts    = addfiles_GetVar(inptr,infiles,"TS")
    ;landfrac = addfiles_GetVar(inptr,infiles,"LANDFRAC")

    ts_land = mask(ts,landfrac.gt.0.5,True)
    sst = mask(ts,landfrac.le.0.33,True)-273

    inptr_gw = addfile(infiles(0),"r")
    gw1  = inptr_gw->gw

    ; calculate restom
    restom = fsnt
    restom = (/fsnt-flnt/)

    ; calculate clearsky restom
    restom_c = fsntc
    restom_c = (/fsntc-flntc/)

    ; weighted global average
    restom_w = wgt_areaave(restom,gw1,1.0,0) 
    restom_w_c = wgt_areaave(restom_c,gw1,1.0,0) 
    flnt_w = wgt_areaave(flnt,gw1,1.0,0) 
    flntc_w = wgt_areaave(flntc,gw1,1.0,0) 
    fsnt_w = wgt_areaave(fsnt,gw1,1.0,0) 
    fsntc_w = wgt_areaave(fsntc,gw1,1.0,0) 
    ts_w = wgt_areaave(ts,gw1,1.0,0) 
    sst_w = wgt_areaave(sst,gw1,1.0,0) 
    ts_land_w = wgt_areaave(ts_land,gw1,1.0,0) 
    if (doTropFlux) then
      flnr_w = wgt_areaave(flnr,gw1,1.0,0) 
      fsnr_w = wgt_areaave(fsnr,gw1,1.0,0) 
    end if
    
    ; create annual average arrays
    restom_ann_w = new(nYears,typeof(restom_w))
    restom_ann_w = 0.

    restom_ann_w_c = new(nYears,typeof(restom_w_c))
    restom_ann_w_c = 0.

    flnt_ann_w = new(nYears,typeof(flnt_w))
    flnt_ann_w = 0.

    fsnt_ann_w = new(nYears,typeof(fsnt_w))
    fsnt_ann_w = 0.

    flntc_ann_w = new(nYears,typeof(flntc_w))
    flntc_ann_w = 0.

    fsntc_ann_w = new(nYears,typeof(fsntc_w))
    fsntc_ann_w = 0.
    
    if (doTropFlux) then
      flnr_ann_w = new(nYears,typeof(flnr_w))
      flnr_ann_w = 0.

      fsnr_ann_w = new(nYears,typeof(fsnr_w))
      fsnr_ann_w = 0.    
    end if

    ts_ann_w = new(nYears,typeof(ts_w))
    ts_ann_w = 0.

    sst_ann_w = new(nYears,typeof(sst_w))
    sst_ann_w = 0.

    ts_land_ann_w = new(nYears,typeof(ts_land_w))
    ts_land_ann_w = 0.

    year = new(nYears,integer)

    if (nYears_a.gt.0) then
      restom_ann_w(0:nYears_a-1) = restom_ann_w_a
      restom_ann_w_c(0:nYears_a-1) = restom_ann_w_c_a
      ts_ann_w(0:nYears_a-1) = ts_ann_w_a
      sst_ann_w(0:nYears_a-1) = sst_ann_w_a
      ts_land_ann_w(0:nYears_a-1) = ts_land_ann_w_a
      year(0:nYears_a-1) = year_a
      flnt_ann_w(0:nYears_a-1) = flnt_ann_w_a
      flntc_ann_w(0:nYears_a-1) = flntc_ann_w_a
      fsnt_ann_w(0:nYears_a-1) = fsnt_ann_w_a
      fsntc_ann_w(0:nYears_a-1) = fsntc_ann_w_a
      if (doTropFlux) then
        flnr_ann_w(0:nYears_a-1) = flnr_ann_w_a
        fsnr_ann_w(0:nYears_a-1) = fsnr_ann_w_a
      end if
    end if

    do iyr=nYears_a,nYears-1
      print((/iyr/))
      iyr2=iyr-nYears_a
      year(iyr)=year_in(iyr2*12)
      do imonth=0,11
        restom_ann_w(iyr) = mon_wgt(imonth)*dim_avg(restom_w(iyr2*12+imonth))+restom_ann_w(iyr)
        print(iyr+" "+imonth+" "+restom_ann_w(iyr)+": "+restom_w(iyr2*12+imonth))
        restom_ann_w_c(iyr) = mon_wgt(imonth)*dim_avg(restom_w_c(iyr2*12+imonth))+restom_ann_w_c(iyr)
        ts_ann_w(iyr) = mon_wgt(imonth)*dim_avg(ts_w(iyr2*12+imonth))+ts_ann_w(iyr)
        sst_ann_w(iyr) = mon_wgt(imonth)*dim_avg(sst_w(iyr2*12+imonth))+sst_ann_w(iyr)
        ts_land_ann_w(iyr) = mon_wgt(imonth)*dim_avg(ts_land_w(iyr2*12+imonth))+ts_land_ann_w(iyr)
        flnt_ann_w(iyr) = mon_wgt(imonth)*dim_avg(flnt_w(iyr2*12+imonth))+flnt_ann_w(iyr)
        flntc_ann_w(iyr) = mon_wgt(imonth)*dim_avg(flntc_w(iyr2*12+imonth))+flntc_ann_w(iyr)
        fsnt_ann_w(iyr) = mon_wgt(imonth)*dim_avg(fsnt_w(iyr2*12+imonth))+fsnt_ann_w(iyr)
        fsntc_ann_w(iyr) = mon_wgt(imonth)*dim_avg(fsntc_w(iyr2*12+imonth))+fsntc_ann_w(iyr)
        if (doTropFlux) then
          flnr_ann_w(iyr) = mon_wgt(imonth)*dim_avg(flnr_w(iyr2*12+imonth))+flnr_ann_w(iyr)
          fsnr_ann_w(iyr) = mon_wgt(imonth)*dim_avg(fsnr_w(iyr2*12+imonth))+fsnr_ann_w(iyr)
        end if
      end do
    end do

    restom_ann_w!0="year"
    restom_ann_w&year=year

    restom_ann_w_c!0="year"
    restom_ann_w_c&year=year

    ts_ann_w!0="year"
    ts_ann_w&year=year

    sst_ann_w!0="year"
    sst_ann_w&year=year

    ts_land_ann_w!0="year"
    ts_land_ann_w&year=year

  ;  year = ispan (year0,nYears+year0-1,1)  

    ; Remove any preexisting file.
    print("  Deleting " + savefile)
    system("/bin/rm " + savefile)

    print(restom_ann_w)

    print("  Creating " + savefile)
    b = addfile(savefile,"c")
    b->year=year
    b->restom_ann_w=restom_ann_w
    b->restom_ann_w_c=restom_ann_w_c
    b->ts_ann_w=ts_ann_w
    b->ts_land_ann_w=ts_land_ann_w
    b->sst_ann_w=sst_ann_w
    b->flnt_ann_w=flnt_ann_w
    b->flntc_ann_w=flntc_ann_w
    b->fsnt_ann_w=fsnt_ann_w
    b->fsntc_ann_w=fsntc_ann_w
    if (doTropFlux) then
      b->flnr_ann_w=flnr_ann_w
      b->fsnr_ann_w=fsnr_ann_w
    end if
    print(b->restom_ann_w)
    delete(b)

  end if

  ;print(restom_ann_w)

  ;print("Average = "+dim_avg(restom_ann_w))


  ; Averages for annual means.

  if (nYears.gt.1) then ; begin plot section
    plot = new(4, graphic)

    resP                       = True
    resP@gsnMaximize           = True            ; expand plot size to maximum
    resP@gsnPaperOrientation   = "landscape"
    resP@txString   = run_name


    res = True
    res@gsnDraw               = False
    res@gsnFrame           = False   ; Don't advance the frame

    res@gsnYRefLine           = 0.              ; reference line   
    res@gsnXYBarChart         = True            ; create bar chart
    res@gsnAboveYRefLineColor = "red"           ; above ref line fill red
    res@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
    res@tiYAxisString = "RESTOM (W/m2)" ; y-axis label      
    res@tiMainString  = run_name ; title


    res1 = True   ; These four lines are redundant.
    res2 = True
    res3 = True
    res4 = True

    res1 = res
    res2 = res
    res3 = res
    res4 = res


    print(year+" "+restom_ann_w+" "+ts_ann_w+" "+ts_land_ann_w+" "+sst_ann_w)

    res1@tiMainString  = "RESTOM" ; title
    res1@gsnYRefLine           = 0.              ; reference line   
    res1@trYMinF = min((/-4,min(restom_ann_w)/))
    res1@trYMaxF = max((/4,max(restom_ann_w)/))
    res1@vpXF              = 0.1    ; X,Y position of first plot
    res1@vpYF              = 0.91
    print("year:"+dimsizes(year))
    print("restom_ann_w:"+dimsizes(restom_ann_w))
    print("RESTOM avg:"+avg(restom_ann_w)+" min:"+min(restom_ann_w)+" max:"+max(restom_ann_w))
    plot(0) = gsn_csm_xy (wks,year,restom_ann_w,res1)

    res2@tiMainString  = "TS" ; title
    res2@tiYAxisString = "TS K" ; y-axis label      
    res2@vpXF              = 0.6    ; X,Y position of second plot
    res2@vpYF              = 0.91
    ;res2@trYMinF = 286.0
    res2@trYMinF = 286.0
    res2@gsnYRefLine           =  res2@trYMinF             ; reference line   
    ;res2@trYMaxF = 289
    res2@trYMaxF = 290.0  ;caf
    print("ts_ann_w:"+dimsizes(ts_ann_w))
    print("TS avg:"+avg(ts_ann_w)+" min:"+min(ts_ann_w)+" max:"+max(ts_ann_w))
    print(res2)
    plot(1) = gsn_csm_xy (wks,year,ts_ann_w,res2)

    res3@tiMainString  = "TS_LAND" ; title
    res3@vpXF              = 0.1    ; X,Y position of third plot
    res3@vpYF              = 0.45
    ;res3@trYMinF = 280.0
    res3@trYMinF = 279.0
    res3@gsnYRefLine           = res3@trYMinF              ; reference line   
    ;res3@trYMaxF = 283
    res3@trYMaxF = 284.0  ;caf
    res3@tiYAxisString = "TS_LAND K" ; y-axis label      
    print("ts_land_ann_w:"+dimsizes(ts_land_ann_w))
    print(min(ts_land_ann_w)+" "+max(ts_land_ann_w))
    plot(2) = gsn_csm_xy (wks,year,ts_land_ann_w,res3)

    res4@tiMainString  = "SST" ; title
    res4@vpXF              = 0.6    ; X,Y position of fourth plot
    res4@vpYF              = 0.45
    res4@tiYAxisString = "SST C" ; y-axis label      
    res4@trYMinF = 16
    res4@trYMaxF = 19
    print("sst_ann_w:"+dimsizes(sst_ann_w))
    print(min(sst_ann_w)+" "+max(sst_ann_w))
    plot(3) = gsn_csm_xy (wks,year,sst_ann_w,res4)

    gsn_panel(wks, plot, (/ 2, 2 /), resP)


    delete(wks)

    print("evince "+pdffile+".pdf")
    evincer=systemfunc("evince " + pdffile + ".pdf &")
    print(evincer)
;  print("mutt mikejmills@mac.com -a "+pdffile+".pdf -s "+pdffile+" < /dev/null")

;  mutter=systemfunc("mutt mikejmills@mac.com -a " + pdffile + ".pdf -s "+pdffile+" < /dev/null")
;  print(mutter)

  end if ; end plot section
end if

