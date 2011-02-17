load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


procedure plotSootOxLoss(historyPattern:string, controlPattern:string)

local ps,p0,hyai,hybi,pm

begin

  ; Some constants we will need:
  pi      = 4.0 * atan(1.0)
  dg2rad  = pi / 180.0   ; rad deg^-1
  g       = 9.806               ; m s^-2 
  re      = 6.37122e06          ; m
  Mair    = 28.964/6.022e23     ; g/molec

  ; Open the set of history files.
  ;
  ; NOTE: We want to work with the entire set of history files, 
  ; so we can do a time evolution.
  print("  Reading data ...")

  historyNames = systemfunc("ls " + historyPath)
  history = addfiles(historyNames, "r")
  ListSetType(history, "cat")
  
  controlNames = systemfunc("ls " + controlPath)
  control = addfiles(controlNames, "r")
  ListSetType(control, "cat")
  
  ; Read in the fields we will need.
  ps = addfiles_GetVar(history, historyNames, "PS")
  p0 = history[0]->P0
  hyai = history[0]->hyai
  hybi = history[0]->hybi
  time = addfiles_GetVar(history, historyNames, "time")
  lev = history[0]->lev
  lat = history[0]->lat
  lon = history[0]->lon
  
  pm = pres_hybrid_ccm(ps,p0,hyai,hybi)
  
  pdel=pm(:,1:dimsizes(lev),:,:)-pm(:,0:dimsizes(lev)-1,:,:)
 
  ; Air density:
  M=addfiles_GetVar(history, historyNames, "M")

  ; Loss rates in molec/cm3/s:
  POx_JO2_x=addfiles_GetVar(history, historyNames, "POx_JO2")
  LOx_CLO_O_x=addfiles_GetVar(history, historyNames, "LOx_CLO_O")
  LOx_HO2_O_x=addfiles_GetVar(history, historyNames, "LOx_HO2_O")
  LOx_HO2_O3_x=addfiles_GetVar(history, historyNames, "LOx_HO2_O3")
  LOx_H_O3_x=addfiles_GetVar(history, historyNames, "LOx_H_O3")
  LOx_NO2_O_x=addfiles_GetVar(history, historyNames, "LOx_NO2_O")
  LOx_O3_O_x=addfiles_GetVar(history, historyNames, "LOx_O3_O")
  LOx_OH_O_x=addfiles_GetVar(history, historyNames, "LOx_OH_O")
  LOx_O_O_x=addfiles_GetVar(history, historyNames, "LOx_O_O")

  POx_JO2_c=addfiles_GetVar(control, controlNames, "POx_JO2")
  LOx_CLO_O_c=addfiles_GetVar(control, controlNames, "LOx_CLO_O")
  LOx_HO2_O_c=addfiles_GetVar(control, controlNames, "LOx_HO2_O")
  LOx_HO2_O3_c=addfiles_GetVar(control, controlNames, "LOx_HO2_O3")
  LOx_H_O3_c=addfiles_GetVar(control, controlNames, "LOx_H_O3")
  LOx_NO2_O_c=addfiles_GetVar(control, controlNames, "LOx_NO2_O")
  LOx_O3_O_c=addfiles_GetVar(control, controlNames, "LOx_O3_O")
  LOx_OH_O_c=addfiles_GetVar(control, controlNames, "LOx_OH_O")
  LOx_O_O_c=addfiles_GetVar(control, controlNames, "LOx_O_O")

  ; Ozone mixing ratios:
  O3_x=addfiles_GetVar(history, historyNames, "O3")
  O3_c=addfiles_GetVar(control, controlNames, "O3")

  



  O3_x_col

  Sum_D_O3      = new((/dimsizes(time),dimsizes(lat)/),float)
  Sum_D_O3!0    = "time"
  Sum_D_O3&time = time
  Sum_D_O3!1    = "lat"
  Sum_D_O3&lat = lat

  Sum_D_O3_O=Sum_D_O3
  Sum_D_NO2_O=Sum_D_O3
  Sum_D_ClO_O=Sum_D_O3
  Sum_D_HO2_O=Sum_D_O3
  Sum_D_HO2_O3=Sum_D_O3
  Sum_D_H_O3=Sum_D_O3
  Sum_D_OH_O=Sum_D_O3
  Sum_D_O_O=Sum_D_O3
  Sum_D_JO2=Sum_D_O3
  Sum_D_Chem=Sum_D_O3

  Sum_O3_O=Sum_D_O3
  Sum_NO2_O=Sum_D_O3
  Sum_ClO_O=Sum_D_O3
  Sum_HO2_O=Sum_D_O3
  Sum_HO2_O3=Sum_D_O3
  Sum_H_O3=Sum_D_O3
  Sum_OH_O=Sum_D_O3
  Sum_O_O=Sum_D_O3
  Sum_JO2=Sum_D_O3
  Sum_Chem=Sum_D_O3
  Sum_O3=Sum_D_O3

  ;Sum_A=Sum_O3_O
  ;Sum_loss=0.0
  ;Sum_D_O3_glob=0.0

  ;integrate over time
  do t=0,dimsizes(time)-1
    print((/time(t)/))
    ; loop over latitudes
    do j=1,dimsizes(lat)-2
    ; loop over altitudes
      do z=0,dimsizes(lev)-1
      ; loop over latitudes
        do l=0,dimsizes(lon)-1
          D_O3=(O3_x(t,z,j,l)-O3_c(t,z,j,l))*M(t,z,j,l)

          D_JO2=POx_JO2_x(t,z,j,l)-POx_JO2_c(t,z,j,l)
          D_O3_O=LOx_O3_O_x(t,z,j,l)-LOx_O3_O_c(t,z,j,l)
          D_NO2_O=LOx_NO2_O_x(t,z,j,l)-LOx_NO2_O_c(t,z,j,l)
          D_ClO_O=LOx_CLO_O_x(t,z,j,l)-LOx_CLO_O_c(t,z,j,l)
          D_HO2_O=LOx_HO2_O_x(t,z,j,l)-LOx_HO2_O_c(t,z,j,l)
          D_HO2_O3=LOx_HO2_O3_x(t,z,j,l)-LOx_HO2_O3_c(t,z,j,l)
          D_H_O3=LOx_H_O3_x(t,z,j,l)-LOx_H_O3_c(t,z,j,l)
          D_OH_O=LOx_OH_O_x(t,z,j,l)-LOx_OH_O_c(t,z,j,l)
          D_O_O=LOx_O_O_x(t,z,j,l)-LOx_O_O_c(t,z,j,l)

          Sum_D_Chem(t,j) = Sum_D_Chem(t,j) + D_JO2 - (D_O3_O + D_NO2_O + D_ClO_O + D_HO2_O + D_HO2_O3 + D_H_O3 + D_OH_O + D_O_O)

    ; eliminate gridpoints where the sign of the change in loss is the same as the sign of the change in O3:      
    ;      if (Sum_loss*D_O3.lt.0) then 

          Sum_D_O3_O(t,j)=Sum_D_O3_O(t,j)+D_O3_O
          Sum_D_NO2_O(t,j)=Sum_D_NO2_O(t,j)+D_NO2_O
          Sum_D_ClO_O(t,j)=Sum_D_ClO_O(t,j)+D_ClO_O
          Sum_D_HO2_O(t,j)=Sum_D_HO2_O(t,j)+D_HO2_O
          Sum_D_HO2_O3(t,j)=Sum_D_HO2_O3(t,j)+D_HO2_O3
          Sum_D_H_O3(t,j)=Sum_D_H_O3(t,j)+D_H_O3
          Sum_D_OH_O(t,j)=Sum_D_OH_O(t,j)+D_OH_O
          Sum_D_O_O(t,j)=Sum_D_O_O(t,j)+D_O_O
          Sum_D_JO2(t,j)=Sum_D_JO2(t,j)+D_JO2

          Sum_D_O3(t,j)=Sum_D_O3(t,j)+D_O3

    ;      end if

        end do
      end do
    end do
  end do

  Sum_O3_O(0,:)=Sum_D_O3_O(0,:)
  Sum_NO2_O(0,:)=Sum_D_NO2_O(0,:)
  Sum_ClO_O(0,:)=Sum_D_ClO_O(0,:)
  Sum_HO2_O(0,:)=Sum_D_HO2_O(0,:)
  Sum_HO2_O3(0,:)=Sum_D_HO2_O3(0,:)
  Sum_H_O3(0,:)=Sum_D_H_O3(0,:)
  Sum_OH_O(0,:)=Sum_D_OH_O(0,:)
  Sum_O_O(0,:)=Sum_D_O_O(0,:)
  Sum_JO2(0,:)=Sum_D_JO2(0,:)
  Sum_Chem(0,:) = Sum_D_Chem(0,:)

  ;integrate over time
  do t=1,dimsizes(time)-1
    print((/time(t)/))
    ; loop over latitudes
    do j=1,dimsizes(lat)-2

          Sum_O3_O(t,j)=Sum_D_O3_O(t,j)+Sum_O3_O(t-1,j)
          Sum_NO2_O(t,j)=Sum_D_NO2_O(t,j)+Sum_NO2_O(t-1,j)
          Sum_ClO_O(t,j)=Sum_D_ClO_O(t,j)+Sum_ClO_O(t-1,j)
          Sum_HO2_O(t,j)=Sum_D_HO2_O(t,j)+Sum_HO2_O(t-1,j)
          Sum_HO2_O3(t,j)=Sum_D_HO2_O3(t,j)+Sum_HO2_O3(t-1,j)
          Sum_H_O3(t,j)=Sum_D_H_O3(t,j)+Sum_H_O3(t-1,j)
          Sum_OH_O(t,j)=Sum_D_OH_O(t,j)+Sum_OH_O(t-1,j)
          Sum_O_O(t,j)=Sum_D_O_O(t,j)+Sum_O_O(t-1,j)
          Sum_JO2(t,j)=Sum_D_JO2(t,j)+Sum_JO2(t-1,j)
          Sum_Chem(t,j) = Sum_D_Chem(t,j) + Sum_Chem(t-1,j) 
    end do
  end do

  fo = addfile("SootOxLoss/2025-05-15-00000.nc","c")

  fo->Sum_D_O3=Sum_D_O3
  fo->Sum_D_Chem=Sum_D_Chem
  fo->Sum_Chem=Sum_Chem

  fo->Sum_D_O3_O=Sum_D_O3_O
  fo->Sum_D_NO2_O=Sum_D_NO2_O
  fo->Sum_D_ClO_O=Sum_D_ClO_O
  fo->Sum_D_HO2_O=Sum_D_HO2_O
  fo->Sum_D_HO2_O3=Sum_D_HO2_O3
  fo->Sum_D_H_O3=Sum_D_H_O3
  fo->Sum_D_OH_O=Sum_D_OH_O
  fo->Sum_D_O_O=Sum_D_O_O
  fo->Sum_D_JO2=Sum_D_JO2

  fo->Sum_O3_O=Sum_O3_O
  fo->Sum_NO2_O=Sum_NO2_O
  fo->Sum_ClO_O=Sum_ClO_O
  fo->Sum_HO2_O=Sum_HO2_O
  fo->Sum_HO2_O3=Sum_HO2_O3
  fo->Sum_H_O3=Sum_H_O3
  fo->Sum_OH_O=Sum_OH_O
  fo->Sum_O_O=Sum_O_O
  fo->Sum_JO2=Sum_JO2

  print(max((/Sum_Chem/)))
  print(max((/Sum_D_O3/)))

  ;************************************************
  ; to plot multiple lines, you must put them into
  ; a mulidimensional array
  ;************************************************
   data      = new((/dimsizes(time),dimsizes(lat)/),float)

   data(0,:) = Sum_Chem(:,{42})
   data(1,:) = Sum_D_O3(:,{42})

  ;************************************************
  ; plotting parameters
  ;************************************************
   wks   = gsn_open_wks ("ps","xy")                  ; open workstation

   res                   = True                      ; plot mods desired
   res@tiMainString      = "Mulitple XY plot"        ; add title

  ; note a similiar resource is xyLineThicknessF, which will effect all
  ; the lines in the array.
   res@xyLineThicknesses = (/1.0,2.0/)               ; make 2nd lines thicker

   res@xyLineColors      = (/"blue","red"/)          ; change line color

   plot  = gsn_csm_xy (wks,time,data,res) ; create plot


;do j=1,dims(2)-2

;  Area=cos(dble2flt(lat(j))*pi/180.)
  
;  print(Area)

;  Sum_O3_O=Sum_O3_O+Sum_D_O3_O(j)*Area

;  Sum_NO2_O=Sum_NO2_O+Sum_D_NO2_O(j)*Area

;  Sum_ClO_O=Sum_ClO_O+Sum_D_ClO_O(j)*Area

;  Sum_HO2_O=Sum_HO2_O+Sum_D_HO2_O(j)*Area

;  Sum_HO2_O3=Sum_HO2_O3+Sum_D_HO2_O3(j)*Area

;  Sum_H_O3=Sum_H_O3+Sum_D_H_O3(j)*Area

;  Sum_OH_O=Sum_OH_O+Sum_D_OH_O(j)*Area

;  Sum_O_O=Sum_O_O+Sum_D_O_O(j)*Area
  
;  Sum_A=Sum_A+Area

;  Sum_D_O3_glob=Sum_D_O3_glob+Sum_D_O3(j)*Area
  
;end do
;Sum_O3_O=Sum_O3_O/Sum_A
;Sum_NO2_O=Sum_NO2_O/Sum_A
;Sum_ClO_O=Sum_ClO_O/Sum_A
;Sum_HO2_O=Sum_HO2_O/Sum_A
;Sum_HO2_O3=Sum_HO2_O3/Sum_A
;Sum_H_O3=Sum_H_O3/Sum_A
;Sum_OH_O=Sum_OH_O/Sum_A
;Sum_O_O=Sum_O_O/Sum_A
;Sum_D_O3_glob=Sum_D_O3_glob/Sum_A

;print(Sum_A)
;print(Sum_O3_O)
;print(Sum_NO2_O)
;print(Sum_ClO_O)
;print(Sum_HO2_O)
;print(Sum_HO2_O3)
;print(Sum_H_O3)
;print(Sum_OH_O)
;print(Sum_O_O)
;print(Sum_D_O3_glob)
