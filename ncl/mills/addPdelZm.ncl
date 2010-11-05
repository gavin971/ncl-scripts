undef("addPdelZm")
procedure addPdelZm(filePath:string)

local i, j, k, b, pdel1, pdel, pdeldry

begin
  
  ; Open the initial conditions file.
  icFile = addfile(filePath, "w")
  
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
 
  ; Calculate the air pressures at the interfaces.
  pdel1 = dpres_hybrid_ccm(PS, P0, hyai ,hybi)
  print("pdel1: "+dimsizes(pdel1))
  pdel1!0="lev"
  pdel1!1="time" 
  pdel1!2="lat"  
  pdel1&lev = lev
  pdel1&lat = lat
  pdel1&time = time
  
  pdel=pdel1(time|:,lev|:,lat|:)
  pdel@long_name="Pressure difference between levels"
  pdel@units="Pa"
  
  pdeldry = pdel ; Get attributes
  pdeldry(:, :, :) = pdeldry(:, :, :)* (1. - Q(time|:, lev|:, lat|:))
  pdeldry@long_name="Dry pressure difference between levels"

  icFile->PDEL=pdel  
  icFile->PDELDRY=pdeldry 
    
end
