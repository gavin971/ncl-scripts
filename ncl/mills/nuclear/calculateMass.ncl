load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("calcPdel")
function calcPdel(historyPath:string)

local ps,p0,hyai,hybi,pm

begin

  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, 
  ; so we can do a time evolution.
  print(" calcPdel: Reading data ...")

  historyNames = systemfunc("ls " + historyPath)
  print("historyNames:"+dimsizes(historyNames))
  history = addfiles(historyPath+historyNames, "r")
  ListSetType(history, "cat")
  
  ; Read in the fields we will need.
  print(" calcPdel: Reading PS")
  ps = addfiles_GetVar(history, historyPath+historyNames, "PS")
  print("ps:"+dimsizes(ps))
  rank=dimsizes(dimsizes(ps))
  p0 = history[0]->P0
  hyai = history[0]->hyai
  hybi = history[0]->hybi
  print(" calcPdel: Reading time")
  time = addfiles_GetVar(history, historyPath+historyNames, "time")
  print("time:"+dimsizes(time))
  lev = history[0]->lev
  lat = history[0]->lat
  if (rank.eq.3) then
    lon = history[0]->lon
  end if
  
  pm = pres_hybrid_ccm(ps,p0,hyai,hybi)
  print("pm:"+dimsizes(pm))
  
;  print(pm(0,:,0,0))
  
  if (rank.eq.3) then
    pdel=pm(:,1:dimsizes(lev),:,:)-pm(:,0:dimsizes(lev)-1,:,:)
  else
    pm!0="lev"
    pm!1="time"
    pm!2="lat"
    pdel=pm(time|:,lev|1:dimsizes(lev),lat|:)-pm(time|:,lev|0:dimsizes(lev)-1,lat|:)
    pdel!0="time"
    pdel!1="lev"
    pdel!2="lat"
  end if
  print("pdel:"+dimsizes(pdel))
  
  print("returning from calcPdel")
;  print(pdel(0,:,0,0))
  
  return(pdel)
  
end



function calculateMass(historyPath:string, tracerNames[*]:string, opt:logical)

; opt: optional argument, true if IC file is to be added via opt@icfile

local pi, dg2grad, g, re, historyNames, hists, lev, lat, lon, nLevs, nLats, nLons, \
  tracerMMR, cellMass, columnMass, dlat, dlon, dx, dy, dydx, wgt, mass, i, j, t, mass

begin

  ; Some constants we will need:
  pi      = 4.0 * atan(1.0)
  dg2rad  = pi / 180.0   ; rad deg^-1
  g       = 9.806               ; m s^-2 
  re      = 6.37122e06          ; m

  nTracers    = dimsizes(tracerNames)
  lastTotal   = 0
  
  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, so we can do
  ; a time evolution of the mass.
  print(" calculateMass: Reading data ...")
  historyNames = systemfunc("ls " + historyPath)
  history = addfiles(historyPath+historyNames, "r")
  ListSetType(history, "cat")
;  histtime=history[:]&time
  histtime=addfiles_GetVar(history, historyPath+historyNames, "time")
;  print(histtime)
  
  ; Read in the fields we will need.
;  pdeldry = addfiles_GetVar(history, historyNames, "PDELDRY")
  print("calling calcPdel")
  pdeldry = calcPdel(historyPath)
  print("returned from calcPdel")
  
  rank=dimsizes(dimsizes(pdeldry))
  print("rank="+rank)

  ; Create the output variable.
  nTimes = dimsizes(histtime)
  print("nTimes="+nTimes)
  print("opt="+opt)
  if (opt) then
    nTimesMass=nTimes+1
    print("nTimesMass="+nTimesMass)
;    print(histtime)
    time=new((/nTimesMass/), typeof(histtime))
    print("time(0)=0")
    time(0)=0.0
    print("time(1:)=histtime")
    time(1:)=histtime
    print("t1=1")
    t1=1
  else
    nTimesMass=nTimes
    time=histtime
    t1=0 
  end if
  
  print("nTimes:"+nTimes+" nTimesMass:"+nTimesMass)
  mass = new((/ nTracers, nTimesMass /), "float")
  mass!0 = "tracer"
  mass&tracer = tracerNames
  mass!1 = "time"
  mass&time = time
  mass@units = "kg"

  do t = 0, nTracers - 1
  
    ; Check for a special tracer name that totals all the other mass.
    print (t+" "+tracerNames(t))
    
    if (tracerNames(t) .eq. "TOTAL") then

      print("  Calculating subtotal ...")
      mass(t, t1:) = mass(lastTotal, :)
      
      do t2 = lastTotal + 1, t-1
        mass(t, t1:) = mass(t, :) + mass(t2, :) 
      end do

    else
    
      print("  Calculating mass for " + tracerNames(t) + " ...")
      tracerMMR = addfiles_GetVar(history, historyPath+historyNames, tracerNames(t))
  
      ; Determine some of the dimensions.
      lev = tracerMMR&lev
      lat = tracerMMR&lat
      nLevs = dimsizes(lev) 
      nLats = dimsizes(lat) 
      if (rank.eq.4) then
        lon = tracerMMR&lon
        nLons = dimsizes(lon)
      end if
    
      ; For each grid point calculate the mass of tracer in the column.
      ;
      ; The amount of TRACER is specifed as a mass mixing ratio (kg/kg) to dry air.
      ; Assuming hydrostatic balance, the column mass at a pressure level is:
      ;    M = P / g
      ; The dry (i.e. no water vapor) pressure drop across the grid point is given
      ; by the field PDELDRY, so the mass of the column is the sum over all the
      ; levels of:
      ;    MASScolumn = SUMlev((PDELDRY / g) * TRACER)
      ;
      ; NOTE: By default, PDELDRY is output to the history files as an average value,
      ; but to check mass conservation, you want to work with instantaneous values.
      print("    cell mass ...")
      print("tracerMMR:"+dimsizes(tracerMMR))
      print("pdeldry:"+dimsizes(pdeldry))
      print("g:"+dimsizes(g))
      cellMass = (pdeldry / g) * tracerMMR
      cellMass!0 = "time"
      cellMass!1 = "lev"
      cellMass&lev = lev
      cellMass!2 = "lat"
      cellMass&lat = lat
      if (rank.eq.4) then      
        cellMass!3 = "lon"
        cellMass&lon = lon
      end if
    
      print("    column mass ...")
      if (rank.eq.4) then      
        columnMass = dim_sum_Wrap(cellMass(time|:, lat|:, lon|:, lev|:))
      else
        columnMass = dim_sum_Wrap(cellMass(time|:, lat|:, lev|:))
      end if
          
      ; Calculate the total mass of the tracer.
      ;
      ; The tracer mass at a surface grid point is the column mass multipled by
      ; the surface area of the grid point:
      ;    MASSgrid = GRIDarea * MASScolumn
      ; The total mass can then be computed adding up the mass of trace from each
      ; grid point:
      ;    MASStotal = SUMlat,lon(MASSgrid)
      
      ; NOTE: Use the wgt_areasum2 function to calculate the sum. The example2 at 
      ;   http://www.ncl.ucar.edu/Document/Functions/Built-in/wgt_areasum2.shtml
      ; shows how to do this.
      
      ; Determine a weighting factor for the sum proportional to the area of the
      ; grid point.
      
      print("    total mass ...")
      if (rank.eq.4) then      
        dlon = fabs(lon(1) - lon(0))
      end if
      dlat = fabs(lat(1) - lat(0))
      
      ; We need to do some special handling to calculate these values to
      ; compensate for the poles and the equator.
      dlat2 = dlat / 2.
      dh = new((/ nLats /), typeof(lat))
    
      do i = 0, nLats - 1
        
        ; The normal handling, +- dlat / 2 in width.
        sdlat = sin((lat(i) + dlat2) * dg2rad) - sin((lat(i) - dlat2) * dg2rad)
        
        ; At the poles, we only want to dlat lat in one direction.
        if (lat(i) .eq. -90.)  then
          sdlat = 1 + sin((lat(i) + dlat2) * dg2rad)
        end if
    
        if (lat(i) .eq. 90.) then
          sdlat =  1 - sin((lat(i) - dlat2) * dg2rad)
        end if
        
        ; At the equator, the peak value is at the equator, so 2 times one side.
        if (lat(i) .eq. 0.)
          sdlat = 2. * (sin((lat(i) + dlat2) * dg2rad))
        end if
    
        dh(i) = re * sdlat
      end do
      
      if (rank.eq.4) then      
        dx = ((2. * pi) / nLons) * re
        dhdx = dh * dx

        wgt = new((/ nLats, nLons /), double)
        wgt = conform(wgt, dhdx, 0)
        mass(t, :) = wgt_areasum2(columnMass, wgt, 0)
      else
        print("columnMass:"+dimsizes(columnMass))
        print(columnMass!0)
        print(columnMass!1)
        wgt = dh * 2. * pi * re
        print("wgt:"+dimsizes(wgt))
        print("wgt sum:"+dim_sum_Wrap(wgt))
        do i = 0, nTimes - 1
          do j = 0, nLats - 1
            columnMass(lat|j,time|i)=doubletofloat(wgt(j))*columnMass(lat|j,time|i)
          end do
        end do
        mass(t, t1:) = (/dim_sum_Wrap(columnMass(time|:,lat|:))/)
      end if
      
    end if
    
    if (opt) then

      icfile=addfile(opt@icfile,"r")
      tracericMMR=icfile->$tracerNames(t)$
      if (rank.eq.4) then
        cellMass3D = (pdeldry(time|0, lat|:, lev|:) / g) * tracericMMR(time|0, lat|:, lev|:)
        cellMass_zm = dim_avg_Wrap(cellMass3D(lat|:, lev|:, lon|:))
      else
        tracer_zm = dim_avg_Wrap(tracericMMR(time|0, lat|:, lev|:, lon|:))
        cellMass_zm = (pdeldry(time|0, lat|:, lev|:) / g) * tracer_zm (lat|:, lev|:)     
      end if
      print("cellMass_zm:"+dimsizes(cellMass_zm))
      cellMass_zm!0 = "lat"
      cellMass_zm!1 = "lev"
      columnMass_zm = dim_sum_Wrap(cellMass_zm(lat|:, lev|:))
      wgt_zm = dh * 2. * pi * re
      massic=dim_sum_Wrap(doubletofloat(wgt(:))*columnMass_zm(:))
      print(massic)
      mass(t,0)=doubletofloat(massic)
      
    end if
    
  end do

;  print(mass&time)
  
  return(mass)
end

