undef("addPdelZm")
procedure addPdelZm(icFile:file,opt:logical)

local i, j, k, b, pdel1, pdel, pdeldry

begin
   notime=False
  ; Open the initial conditions file.
;  icFile = addfile(filePath, "w")
  
  ; Read in the dimensions.
  lev 	= icFile->lev
  lat 	= icFile->lat
  time 	= icFile->time
  
  ; Do we need to iterate (or do array math) to get ptdeldry at each
  ; grid location and vertical level.
  hyai 	= icFile->hyai
  hybi 	= icFile->hybi
  P0 	= icFile->P0
  PS 	= icFile->PS
  Q     = icFile->Q
  
  if (dimsizes(dimsizes(PS)) .eq. 1) then
    notime=True
    PS2=new((/1,dimsizes(PS)/),float)
    PS2(0,:)=PS
    PS2!0="time"
    delete(PS)
    PS=PS2
  end if
 
  ; Calculate the air pressures at the interfaces.
  pdel1 = dpres_hybrid_ccm(PS, P0, hyai ,hybi)
;  print("pdel1: "+dimsizes(pdel1))
;  print("lev: "+dimsizes(lev))
  pdel1!0="lev"
  pdel1!1="time" 
  pdel1!2="lat"  
  pdel1&lev = lev
  pdel1&lat = lat
  pdel1&time = time
  
  if (notime) then
    pdel=pdel1(time|0,lev|:,lat|:)   
  else
    pdel=pdel1(time|:,lev|:,lat|:)
  end if
  pdel@long_name="Pressure difference between levels"
  pdel@units="Pa"
  
  pdeldry = pdel ; Get attributes
  pdeldry = pdeldry* (1. - Q)
  pdeldry@long_name="Dry pressure difference between levels"

  print("  writing PDEL")
  icFile->PDEL=pdel  
  print("  writing PDELDRY")
  icFile->PDELDRY=pdeldry
  
  if (opt) then
    opt@pdel=pdel
    opt@pdeldry=pdeldry
    opt@notime=notime
  end if
    
end
