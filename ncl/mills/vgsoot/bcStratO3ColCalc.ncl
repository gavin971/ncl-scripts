load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "../addPdelZm.ncl"

undef("calcZdelZavg")
function calcZdelZavg(a:file)

local z3,zdel
begin
 
  print("Reading geopotential heights...")
  z3 = a->Z3
;  print("Z3:"+dimsizes(z3))
  nlev=dimsizes(z3&lev)
  
  print("Calculating zdel...")
  zdel = z3
  
  rank=dimsizes(dimsizes(zdel))
  if (rank.eq.3) then
    zdel(:,0:nlev-2,:)=z3(:,0:nlev-2,:)-z3(:,1:nlev-1,:)
  else
    zdel(0:nlev-2,:)=z3(0:nlev-2,:)-z3(1:nlev-1,:)
  end if
    
  delete(z3)
  delete(nlev)
  
  return(zdel)
    
end

undef("calcStratColZavg_vmr")
function calcStratColZavg_vmr(pdeldry:float, vmr:float)

local g,Mair,partPress,column,vmr,columnAv
begin
  ; Some constants we will need:
  g       = 980.6               ; cm s^-2 
  Mair    = 28.964/6.022e23     ; g/molec
  pi=acos(0.0)*2.0

  ; dP=rho*g*dh=[M]*Mair*g*dh
  ; [M]*dh=dP/(Mair*g)
  partPress=pdeldry*vmr*10./(Mair*g)  ; pdeldry(Pa)*10(dynes/cm2/Pa)
  partPress!0 = "time"
  partPress&time = vmr&time
  partPress!1 = "lev"
  partPress&lev = vmr&lev
  partPress!2 = "lat"
  partPress&lat = vmr&lat
    
  print("  Calculating column from vmr...")
  column = dim_sum_Wrap(partPress(time|:, lat|:, {lev|1:100}))

  return(column)
end

undef("saveZavgCol")
procedure saveZavgCol(historyRun:string, rootPath:string, season:string)

local historyPattern,historyNames,history,T,Tmean,pi,pm,press,temp,tmean

begin

  pi=acos(0.0)*2.0
  
  historyPattern = rootPath+"/"+historyRun+"/h0zm."+season+"/"+historyRun+".cam2.h0zm.20*."+season+".nc"
  historyNames = systemfunc("ls " + historyPattern)
  zmNames = historyNames
  nFiles=dimsizes(historyNames)
   
  i = 0
  do while ((i.lt.nFiles).and.(isfilepresent(zmNames(i))))
    print("File "+i+": "+zmNames(i))
    a = addfile(zmNames(i),"w")
    lat = a->lat
    nlat = dimsizes(lat)
    
    ; pdeldry:
    if (isfilevar(a,"PDELDRY").and.isfilevar(a,"PDEL")) then
      print("  reading PDELDRY...")
      pdeldry=a->PDELDRY
    else
      print("  computing PDELDRY...")
      opt=True
      addPdelZm(a,opt)
      pdeldry=opt@pdeldry
    end if
        
    ; Ozone:
    if (isfilevar(a,"O3stcol")) then
      print("  O3stcol exists, skipping...")
    else
      print("  reading ozone...")
      O3=a->O3

      O3_col=new((/1,nlat/),float)
      O3_col!0 = "time"
      O3_col&time = O3&time
      O3_col!1 = "lat"
      O3_col&lat = O3&lat

      O3_col=calcStratColZavg_vmr(pdeldry, O3)
      O3_col=O3_col/2.687e16 ; Convert to Dobson Units
      O3_col@units="Dobson Units"
      O3_col@long_name="Stratospheric O3 Column, 100hPa>p>1hPa"

      print("  saving ozone column...")
      a->O3stcol=O3_col
      delete(O3_col)
    end if

    i = i + 1
  end do
  
end


;rootPath="/data5/mmills/wa319/carma-soot"
rootPath="/Volumes/Data/Models/waccm/run/carma-soot"

historyRun="bcrckt03"
;historyRun="bccntl02"
;season="DJF"
season="JJA"

saveZavgCol(historyRun, rootPath, season)

