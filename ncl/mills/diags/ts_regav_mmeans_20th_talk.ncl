; Prior to running, do:
;   ncrcat -v TS <h0.files*.nc> <TS_timeseries.nc>
; and have that as your addfile input in ncl.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

host=getenv("HOSTNAME")

if (host.eq."cmip5") then
  in_dir = "/data/mmills/"
else
  in_dir = "/data3/mmills/"
end if

;year0 = (/1850, 1960, 1955, 1850, 1850, 1850/)
;yearmax=(/2005, 2049, 2031, 2005, 2059, 2005/)
;run_names = (/"b40.20th.track1.2deg.wcm.007","b40.1955-2005.2deg.wcm.002",\
;              "b40.1955-2005.2deg.wcm.003",\
;              "b40.20th.2deg.wset.001","b40.rcp4_5.1deg.ave","Obs"/)
;case_names = (/"WACCM4 20th Century","WACCM4 Ensemble 1","WACCM4 Ensemble 2",\
;               "CAM4 2deg","CAM4 1deg 5-Ensemble Average","HadCRU"/)
;xyLineColors      = (/"blue","darkgreen","darkgreen","red","red","black"/)    
;xyDashPatterns = (/0,0,1,0,1,1/)

;year0 = (/1850, 1850, 1850, 1960, 1955, 1955, 1850/)
;yearmax=(/2005, 1957, 1942, 2072, 2051, 2005, 2005/)
;run_names = (/"b40.20th.track1.2deg.wcm.007",\
;              "br04.20th.track1.2deg.wcm.pd.001",\
;              "br04.20th.track1.2deg.wcm.pd.002",\
;              "b40.1955-2005.2deg.wcm.002",\
;              "b40.1955-2005.2deg.wcm.003",\
;              "b40.1955-2005.2deg.wcm.004",\
;              "Obs"/)
;case_names = (/"WACCM4 20th Century 1",\
;               "WACCM4 20th Century 2",\
;               "WACCM4 20th Century 3",\
;               "WACCM4 CMIP5 1",\
;               "WACCM4 CMIP5 2",\
;               "WACCM4 CMIP5 3",\
;               "HadCRU"/)
;xyLineColors      = (/"blue","turquoise","orange","red","darkgreen","chartreuse","black"/)    
;xyDashPatterns = (/0,0,0,0,0,0,1/)

year0 = (/1850, 1960, 2005, 2005, 1850/);, 1850/)
yearmax=(/2005, 2099, 2099, 2099, 2059/);, 2005/)
run_names = (/"b40.20th.track1.2deg.wcm.avg",\
              "b40.1955-2005.2deg.wcm.avg",\
              "b40.rcp8_5.2deg.wcm.003",\
              "b40.rcp2_6.2deg.wcm.002",\
              "b40.rcp4_5.1deg.ave"/);,\
;              "Obs"/)
case_names = (/"WACCM4 20thC Ave(3)",\
               "WACCM4 RCP4.5 Ave(3)",\
               "WACCM4 RCP8.5",\
               "WACCM4 RCP2.6",\
               "CAM4 1deg RCP4.5 Ens Ave(5)"/);,\
;               "HadCRU"/)
xyLineColors      = (/"blue","darkgreen","red","black","darkgreen"/)    
;xyDashPatterns = (/0,0,0,1/)
xyDashPatterns = (/0,0,0,0,1/)

con_names = (/"b40.1850.track1.1deg.006","Obs"/)

maxyear=max(yearmax)
minyear=min(year0)
nyearsmax = maxyear-minyear+1
filt121 = True
nyrunave = 11
omask0 = False
;omask0 = True
;norm_plot = True
norm_plot = False
norm_plot_yrs_abs = (/1961,1990/)

;var_names = (/"LWCF","RESTOM","FSNT","FLNT"/)
var_names = (/"TS"/)
nyears = yearmax - year0 + 1

con_tav = (/999.99,287.657,287.206,0./)

ptype = (/"Global","Land","Northern Hemisphere","Southern Hemisphere"/)

lrmcon = False

region = "Global"
mask_obs = False
ldtrend = False

ps_pref = "LBNL"

xcor = (/0.12,0.63,0.12,0.63/)
ycor = (/0.70,0.70,0.31,0.31/)

yr_weight = new(nyrunave,"float")
yr_weight = 1.

nplots = 4

if (lrmcon) then
  ps_pref = ps_pref+"_rmcon"
end if
;;;;;;;;

;obs_file = "/project/amp/rneale/data/HADCRU/anna_obs_abs.nc"
obs_file = "$DATA_PATH/HADCRU/HadCRUT3v.nc"

;;;;;;;;;;;;;;;;;;;;

ncases = dimsizes(case_names)
nvarnames = dimsizes(var_names)
year = ispan (minyear,maxyear,1)

mon_wgt = (/0.08493150770664215,0.07671232521533966,0.08493150770664215,0.08219178020954132,0.08493150770664215,0.08219178020954132,0.08493150770664215,0.08493150770664215,0.08219178020954132,0.08493150770664215,0.08219178020954132,0.08493150770664215/)

Lf = 3.337e5
Lv = 2.501e6    

; Calculate normalizing years.

norm_plot_yrs = norm_plot_yrs_abs - 1850


if (region.eq."Global") then
  lats = -90.
  latn = 90.
  lonw = 0.
  lone = 360.
  reg_title = "Global ("+lats+" - "+latn+"N, "+lonw+" - "+lone+"E)" 
end if
if (region.eq."NHem") then
  lats = 60.
  latn = 90.
  lonw = 0.
  lone = 360.
  reg_title = "Northern of 60N ("+lats+" - "+latn+"N, "+lonw+" - "+lone+"E)" 
  con_tav = (/287.657,287.694,287.206,0./)
end if
if (region.eq."SHem") then
  lats = -90.
  latn = -60.
  lonw = 0.
  lone = 360.
  reg_title = "South of 60S ("+lats+" - "+latn+"N, "+lonw+" - "+lone+"E)" 
  con_tav = (/287.657,287.694,287.206,0./)
end if

txres = True
txres@txFontHeightF = 0.012

pan = True
pan@gsnMaximize = True
pan@gsnPaperOrientation = "landscape"

res = True
res@gsnDraw = False
res@gsnFrame = False

;res@trYMinF = -2.5
;res@trYMaxF = 2.5
;res@trXMaxF = maxyear+1
;res@trXMaxF = 2050
res@trXMaxF = 2099
res@gsnYRefLine           = 0.              ; reference line   
res@gsnYRefLineDashPattern = 1.
res@xyLineColors   = xyLineColors
res@xyDashPatterns = xyDashPatterns
res@xyLineThicknessF = 2
 

res@pmLegendDisplayMode    = "Never"            ; turn on legend
res@lgBoxLinesOn = False
res@pmLegendSide           = "Top"               ; Change location of 
res@pmLegendParallelPosF   = .25                 ; move units right
res@pmLegendOrthogonalPosF = -0.3                ; move units down

res@pmLegendWidthF         = 0.15                ; Change width and
res@pmLegendHeightF        = 0.15                ; height of legend.
res@lgLabelFontHeightF     = .015                 ; change font height

;res@vpYF            = 0.3
;res@vpHeightF       = .2                    ;
res@vpWidthF        = .9
res@xyExplicitLegendLabels = case_names
res@lgTitleFontHeightF = 0.015
res@tiXAxisString = "Year"
;res@tiYAxisString = "Anomaly"
;

do ivar=0,nvarnames-1

  plot1 = new(nplots,graphic)
  plot2 = new(nplots,graphic)
  
  var_name = var_names(ivar)
  print(var_name+"")
  if (any(var_name.eq.(/"LHFLX","FLNT","FSNT","FSNTC","FLNTC","SWCF","LWCF","RESTOM"/))) then
    res@tiYAxisString = "Flux anomaly (Wm~S~-2~N~)"
    res@trYMinF = -2.5
    res@trYMaxF = 2.5
  end if

  if (any(var_name.eq.(/"FSNS","FSDS","FSDSC","FLNS","FLDS"/))) then
    res@tiYAxisString = "Flux anomaly (Wm~S~-2~N~)"
    res@trYMinF = -5.
    res@trYMaxF = 5
  end if

  if (any(var_name.eq.(/"TS","TREFHT"/))) then
    if(norm_plot) then
      res@tiYAxisString = "TS anomaly (K)"
      res@trYMinF = -1.2
      res@trYMaxF = 1.8
    else
      res@tiYAxisString = "TS (K)"
    end if
  end if
  if (any(var_name.eq.(/"AODVIS"/))) then
    res@tiYAxisString = "AOD anomaly"
    res@trYMinF = -0.01
    res@trYMaxF = 0.05
  end if
  if (any(var_name.eq.(/"PRECT"/))) then
    res@tiYAxisString = "Prect. anomaly (mm/day)"
    res@trYMinF = -0.3
    res@trYMaxF = 0.3
  end if
  if (any(var_name.eq.(/"TMQ"/))) then
    res@tiYAxisString = "PW anomaly (mm)"
    res@trYMinF = -1.5
    res@trYMaxF = 1.5
  end if
  if (any(var_name.eq.(/"CLDLOW","CLDHGH","CLDMED","CLDTOT"/))) then
    res@tiYAxisString = "Cloud anomaly (fraction)"
    res@trYMinF = -0.015
    res@trYMaxF = 0.015
  end if


  pan@txString = var_name

  var_ann_w = new((/ncases,nyearsmax/),"float")
  var_c_w  = new(ncases,"float")

  filename=  "$IMG_PATH/"+ps_pref+"_20th_gmean_ts_"+var_name
  print(filename)
  wks   = gsn_open_wks ("ps", filename)
  wks_s = gsn_open_wks ("ps", filename+"_smooth")

  pi_ave = new(ncases,"string")
  var_ann_w = 0.

  do ip=0,nplots-1 ; 4 plots (global, land, nhem, shem)
; Read and weight
    print(ptype(ip)+"")
    if (ip.eq.0) then ; Global and land
      lats = -90.
      latn = 90.
      lonw = 0.
      lone = 360.
    end if
    if (ip.eq.2) then ; NHem
      lats = 0.
      latn = 90.
      lonw = 0.
      lone = 360.
    end if
    if (ip.eq.3) then ; S. Hemisphere
      lats = -90
      latn = 0.
      lonw = 0.
      lone = 360.
    end if


    do icase=0,ncases-1
      print(run_names(icase)+"")
      
      if (run_names(icase).ne."Obs")
        
        if (var_name.eq."RESTOM") then
          infiles1 = systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*FSNT.20th.nc")
          infiles2 = systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*FLNT.20th.nc")
          inptr1 = addfile(infiles1,"r")
          inptr2 = addfile(infiles2,"r")
          fsnt    = inptr1->FSNT(:,{lats:latn},{lonw:lone})
          flnt    = inptr2->FLNT(:,{lats:latn},{lonw:lone})
          var_in =  fsnt
          var_in = (/fsnt-flnt/)
          var_in@long_name = "RESTOM"
          delete(fsnt)
          delete(flnt)        
        end if
        if (var_name.eq."PRECT") then
          infiles1 = systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*PRECC.20th.nc")
          infiles2 = systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*PRECL.20th.nc")
          inptr1 = addfile(infiles1,"r")
          inptr2 = addfile(infiles2,"r")
          precc    = inptr1->PRECC(:,{lats:latn},{lonw:lone})
          precl    = inptr2->PRECL(:,{lats:latn},{lonw:lone})
          var_in =  precc
          var_in = (/86400.*1000.*(precc+precl)/)
          var_in@long_name = "Total Precipitation"
          delete(precc)
          delete(precl)  
        end if
        if (.not.isdefined("var_in")) then
          if (var_name.eq."FSDSC".and.icase.eq.1) then
            infiles = systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*FSDS.20th.nc")
            inptr = addfile(infiles,"r")
            var_in = inptr->FSDS(:,{lats:latn},{lonw:lone})
            var_in = var_ann_w@_FillValue
          else
            if (var_name.eq."AODVIS".and.icase.eq.2) then
              infiles = systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*AEROD_v.20th.nc")
              inptr = addfile(infiles,"r")
              var_in = inptr->AEROD_v(:,{lats:latn},{lonw:lone})
            else
              infiles = systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*"+var_name+".20th.nc")
              inptr = addfile(infiles,"r")
              var_in = inptr->$var_name$(:,{lats:latn},{lonw:lone})
            end if
          end if
        end if

        if (var_name.eq."TS".and.lrmcon) then
          infile_c =  systemfunc("ls -1 "+in_dir+con_names(icase)+"/proc/"+con_names(icase)+"_ANN_climo.nc")
          inptr_c = addfile(infile_c,"r")
          var_c = inptr_c->TS(0,{lats:latn},:)
        end if

; Gaussian weights
        infile_gw =  systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*gw*nc")
        inptr_gw = addfile(infile_gw,"r")
        gw  = inptr_gw->gw({lats:latn})

; Land frac
        infile_frac =  systemfunc("ls -1 "+in_dir+run_names(icase)+"/proc/*LANDFRAC*nc")
        inptr_frac = addfile(infile_frac,"r")
        landfrac  = inptr_frac->LANDFRAC(0,{lats:latn},:)
        
        if (var_name.eq."AODVIS") then
          delete(var_in@_FillValue)
          var_in@_FillValue = 1.e36
          if (var_name.eq."AODVIS") then
            res@trYMinF = -0.02
            res@trYMaxF = 0.1
          end if
        end if
; Maskk with obs?
        
        if (mask_obs) then
          
          inptr_o = addfile(obs_file,"r")
          var_ino = inptr_o->temp(:,0,{lats:latn},:)
          lato = var_ino&latitude
          lono = var_ino&longitude
          gwo = latRegWgt(lato,"float",0) 
          
          var_in_ogrid = var_ino
          var_in_ogrid = (/var_ino*0./)
          
;        obsvar = (/mask(obsvar, sst_obs_int.eq.sst_obs_int@_FillValue, False)/)
          im = 0
          print(var_in_ogrid(0,29,:))
          do iyr=0,0
            do imonth=0,1
              var_in_ogrid_tmp = linint2(var_in&lon,var_in&lat,var_in(im,:,:),True,lono,lato,0)
              var_in_ogrid(im,:,:) = var_in_ogrid_tmp
              delete(var_in_ogrid_tmp)
;            var_in_ogrid(im,:,:) = var_in_ogrid(im,:,:)
              var_in_ogrid(im,:,:) =  mask(var_in_ogrid(im,:,:),var_ino(im,:,:).eq.var_ino@_FillValue,False) 
              im = im+1
            end do
          end do
          print("after "+var_in_ogrid(0,29,:))
          
          delete(var_in)
          delete(gw)
          var_in =  var_in_ogrid 
          gw = gwo
          delete(var_in_ogrid)
          delete(gwo)
          
        end if
        
;Weighting coupled runs.
; Land maskk

        if (ip.eq.1) then ; Land
          var_in = mask(var_in,landfrac.gt.0.75,True)
          if ((lrmcon.and.var_name.eq."TS")) then
            var_c  = mask(var_c,landfrac.gt.0.75,True)
          end if
        end if
        delete(landfrac)
        


; Weight timeseries
        
        var_w = wgt_areaave(var_in,gw,1.0,0) 
        if ((lrmcon.and.var_name.eq."TS")) then
          var_c_w(icase) = wgt_areaave(var_c,gw,1.0,0) 
          delete(var_c)
        end if
        
        do iyr=0,nyears(icase)-1
          iyr2=iyr+year0(icase)-1850
          do imonth=0,11
;            print(run_names(icase)+" "+nyears(icase)+" "+year0(icase)+" "+iyr+" "+imonth)
            var_ann_w(icase,iyr2) = mon_wgt(imonth)*dim_avg(var_w(iyr*12+imonth))+var_ann_w(icase,iyr2)
          end do
;          print(run_names(icase)+" iyr:"+iyr+" iyr2:"+iyr2+" var:"+var_ann_w(icase,iyr2))
        end do
        if (ldtrend) then
          var_dtrend = dtrend(var_ann_w(icase,0:49),True)
          printVarSummary(var_dtrend)
          do iyr=0,nyears(icase)-1
            iyr2=iyr+year0(icase)-1850
            var_ann_w(icase,iyr2) = var_ann_w(icase,iyr2)-iyr*var_dtrend@slope
          end do
          delete(var_dtrend)
        end if
        delete(var_in)
        delete(gw)
        delete(var_w)
      else
      
        inptr_obs = addfile(obs_file,"r")
;      var_w = inptr_o->ABS_OBS
        var_in = inptr_obs->temp(:,0,{lats:latn},:)
;      printVarSummary(var_in)
        latobs = var_in&latitude
        gwobs = latRegWgt(latobs,"float",0) 
        
        var_w = wgt_areaave(var_in,gwobs,1.0,0) 
        delete(var_in)
        
      
        do iyr=0,nyears(icase)-1
          iyr2=iyr+year0(icase)-1850
          do imonth=0,11
            var_ann_w(icase,iyr2) = mon_wgt(imonth)*dim_avg(var_w(iyr*12+imonth))+var_ann_w(icase,iyr2)
          end do
        end do
        
        
        if (var_name.ne."TS") then
          var_ann_w(icase,:) = var_ann_w@_FillValue
        end if
        delete(var_w)
        delete(latobs)
        delete(gwobs)
        
      end if
      
;      var_ann_w(icase,nyears(icase)-1:nyearsmax-1) = var_ann_w@_FillValue
      var_ann_w(icase,yearmax(icase)-1850:) = var_ann_w@_FillValue
      var_ann_w(icase,:year0(icase)-1850) = var_ann_w@_FillValue
          
    end do
    
    var_ann_w0=var_ann_w(0,norm_plot_yrs(0):norm_plot_yrs(1))
    
    TS_0=dim_avg(var_ann_w0)
    print("TS_0="+TS_0)
;    case_names(0)=case_names(0)+": TS~B~0~N~="+TS_0
    
; Averages for annual means.

    if (norm_plot) then
      do icase=0,ncases-1
  ;      print(run_names(icase)+" "+norm_plot_yrs(0)+" "+norm_plot_yrs(1))
        print(dimsizes(var_ann_w))
        pi_ave = " (PI ave = "+flt2string(dim_avg(var_ann_w(icase,norm_plot_yrs(0):norm_plot_yrs(1))))+")"
        res@xyExplicitLegendLabels(icase) = case_names(icase)+pi_ave(icase)
        if (lrmcon.and.icase.ne.ncases-1.and.var_name.eq."TS") then
          var_ann_w(icase,:) =  var_ann_w(icase,:)-var_c_w(icase)
        else
;          if (run_names(icase).eq."Obs") then
            var_ann_w(icase,:) =  var_ann_w(icase,:)-dim_avg(var_ann_w(icase,norm_plot_yrs(0):norm_plot_yrs(1)))
;          else
;            var_ann_w(icase,:) =  var_ann_w(icase,:)-TS_0
;          end if
        end if
      end do
    end if
    
    if (var_name.eq."FLNT".or.var_name.eq."FLNTC") then
      var_ann_w = -var_ann_w
    end if

    res@tiMainString = ptype(ip)
;    res@gsnYRefLines = (/0.0,0.0,0.0,0.0/)
;    print(res@xyLineColors)
;    res@gsnAboveYRefLineColors = res@xyLineColors 
;    res@gsnBelowYRefLineColors = res@xyLineColors 
    if (ip.ne.1) then
      plot1(ip) = gsn_csm_xy (wks,int2flt(year),var_ann_w,res)
    ; Add text string to small map.
      txres2               = True
      txres2@txFontHeightF = 0.04
      txres2@txFontColor = "black"
;      txres2@txAngleF      = -45.     ; Rotate text 45 degrees
      textstr="TS~B~0~N~="+TS_0
      text_anno = gsn_add_text(wks,plot1(ip),textstr,1890,1.0,txres2)
    end if

    if (var_name.eq."TS") then
      cor1 = sprintf("%8.3f",escorc(var_ann_w(0,0:nyears(0)-1),var_ann_w(ncases-1,0:nyears(0)-1)))
      cor2 = sprintf("%8.3f",escorc(var_ann_w(1,0:nyears(0)-1),var_ann_w(ncases-1,0:nyears(0)-1)))
      txres@txFontColor = "blue"
;      gsn_text_ndc(wks,cor1,xcor(ip),ycor(ip),txres) 
      txres@txFontColor = "darkgreen"
;      gsn_text_ndc(wks,cor2,xcor(ip),ycor(ip)+0.02,txres) 
    end if
    
    do icase=0,ncases-1
      if (filt121) then
;        do ifpass = 0,nyears(icase)-1
          var_ann_w(icase,:) = wgt_runave(var_ann_w(icase,:),yr_weight, 0)
;        end do
      end if
      
    end do
    
    
    i1 = 2004-1850
    i2 = 1850-1850
    if (ip.ne.1) then
      plot2(ip) = gsn_csm_xy (wks_s,int2flt(year),var_ann_w,res)
    end if
    if (var_name.eq."TS") then
      cor1 = sprintf("%8.3f",escorc(var_ann_w(0,0:nyears(0)-1),var_ann_w(ncases-1,0:nyears(0)-1)))
      cor2 = sprintf("%8.3f",escorc(var_ann_w(1,0:nyears(0)-1),var_ann_w(ncases-1,0:nyears(0)-1)))
      txres@txFontColor = "blue"
;      gsn_text_ndc(wks_s,cor1,xcor(ip),ycor(ip),txres)
      txres@txFontColor = "darkgreen"
;      gsn_text_ndc(wks_s,cor2,xcor(ip),ycor(ip)+0.02,txres) 
end if
    
    var_ann_w = (/0./)

  end do ; Global, land ,nhem, shem

  pan@gsnFrame       = False
  gsn_panel(wks,plot1,(/2,2/),pan)
  gsn_panel(wks_s,plot2,(/2,2/),pan)

   legend = create "Legend" legendClass wks 
    "vpXF"                     : 0.6                    ; orientation on page
    "vpYF"                     : 0.73
    "vpWidthF"                 : 0.08                     ; width
    "vpHeightF"                : 0.1                     ; height
    "lgPerimOn"                : True                   ; no perimeter
    "lgItemCount"              : ncases                      ; how many
    "lgLabelStrings"           : case_names                   ; Labels
    "lgLabelsOn"               : True                   ; no default lables
    "lgMonoDashIndex"          : False
    "lgDashIndexes"            : res@xyDashPatterns
    "lgLineThicknessF"         : res@xyLineThicknessF 
    "lgLabelFontHeightF"       : 0.010                   ; font height
    "lgLineColors"             : res@xyLineColors
    "lgMonoLineLabelFontColor" : False                    ; one label color
    "lgItemOrder"              : ispan(ncases-1,0,1)    ; Reverse order 
    "lgAutoManage"             : False
    "lgLabelJust"              : "CenterLeft"
    end create
   
     legend_s = create "Legend" legendClass wks_s 
    "vpXF"                     : 0.6                    ; orientation on page
    "vpYF"                     : 0.73
    "vpWidthF"                 : 0.08                     ; width
    "vpHeightF"                : 0.1                     ; height
    "lgPerimOn"                : True                   ; no perimeter
    "lgItemCount"              : ncases                      ; how many
    "lgLabelStrings"           : case_names                   ; Labels
    "lgLabelsOn"               : True                   ; no default lables
    "lgMonoDashIndex"          : False
    "lgDashIndexes"            : res@xyDashPatterns 
    "lgLineThicknessF"         : res@xyLineThicknessF
    "lgLabelFontHeightF"       : 0.010                   ; font height
    "lgLineColors"             : res@xyLineColors
    "lgMonoLineLabelFontColor" : False                    ; one label color
    "lgItemOrder"              : ispan(ncases-1,0,1)    ; Reverse order 
    "lgAutoManage"             : False
    "lgLabelJust"              : "CenterLeft"
    end create

    draw(legend)
    draw(legend_s)
     
   
    pan@txString = var_name  
   
    
    frame(wks)
    pan@txString = var_name+" - Smoothed"  
    frame(wks_s)

  delete(plot1)
  delete(plot2)
  delete(var_ann_w)
;  if (ncases.eq.3) then
;    ncases = ncases-1
;  end if
  
end do

;print("TS_0="+TS_0)

print("ps2pdf " + filename + ".ps" + " " + filename + ".pdf")
pstest=systemfunc("ps2pdf " + filename + ".ps" + " " + filename + ".pdf")
print(pstest)

print("mutt -a " + filename + ".pdf mikejmills@mac.com < /dev/null")
mutter=systemfunc("mutt -a " + filename + ".pdf mikejmills@mac.com < /dev/null")
print(mutter)

end
