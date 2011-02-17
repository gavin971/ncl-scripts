load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

x=addfile("OxRateRad80mw/OxRateRad80mw.cam2.h0.2025-05-15-00000.nc","r")
control=addfile("OxRateTest2/OxRateTest2.cam2.h0.2025-05-15-00000.nc","r")
pi=acos(0.0)*2.0

;(time, lev, lat, lon)
time=x->time
lev=x->lev
lat=x->lat
lon=x->lon
t=24
print(time(t))

M=x->M

LOx_CLO_O_x=x->LOx_CLO_O
LOx_HO2_O_x=x->LOx_HO2_O
LOx_HO2_O3_x=x->LOx_HO2_O3
LOx_H_O3_x=x->LOx_H_O3
LOx_NO2_O_x=x->LOx_NO2_O
LOx_O3_O_x=x->LOx_O3_O
LOx_OH_O_x=x->LOx_OH_O
LOx_O_O_x=x->LOx_O_O

LOx_CLO_O_c=control->LOx_CLO_O
LOx_HO2_O_c=control->LOx_HO2_O
LOx_HO2_O3_c=control->LOx_HO2_O3
LOx_H_O3_c=control->LOx_H_O3
LOx_NO2_O_c=control->LOx_NO2_O
LOx_O3_O_c=control->LOx_O3_O
LOx_OH_O_c=control->LOx_OH_O
LOx_O_O_c=control->LOx_O_O

O3_x=x->O3
O3_c=control->O3

dims=dimsizes(LOx_O3_O_x)
Sum_D_O3=O3_x(0,0,:,0)
Sum_D_O3=0.0
Sum_D_O3_O=Sum_D_O3
Sum_D_NO2_O=Sum_D_O3
Sum_D_ClO_O=Sum_D_O3
Sum_D_HO2_O=Sum_D_O3
Sum_D_HO2_O3=Sum_D_O3
Sum_D_H_O3=Sum_D_O3
Sum_D_OH_O=Sum_D_O3
Sum_D_O_O=Sum_D_O3

Sum_O3_O=0.0
Sum_NO2_O=0.0
Sum_ClO_O=0.0
Sum_HO2_O=0.0
Sum_HO2_O3=0.0
Sum_H_O3=0.0
Sum_OH_O=0.0
Sum_O_O=0.0

Sum_A=Sum_O3_O
Sum_loss=0.0
Sum_D_O3_glob=0.0

; loop over latitudes
do j=0,dims(2)-1
; loop over altitudes
  do z=0,dims(1)-1
  ; loop over latitudes
    do l=0,dims(3)-1
      Sum_D_O3(j)=Sum_D_O3(j)+(O3_x(t,z,j,l)-O3_c(t,z,j,l))*M(t,z,j,l)
      Sum_D_O3_O(j)=Sum_D_O3_O(j)+LOx_O3_O_x(t,z,j,l)-LOx_O3_O_c(t,z,j,l)
      Sum_D_NO2_O(j)=Sum_D_NO2_O(j)+LOx_NO2_O_x(t,z,j,l)-LOx_NO2_O_c(t,z,j,l)
      Sum_D_ClO_O(j)=Sum_D_ClO_O(j)+LOx_CLO_O_x(t,z,j,l)-LOx_CLO_O_c(t,z,j,l)
      Sum_D_HO2_O(j)=Sum_D_HO2_O(j)+LOx_HO2_O_x(t,z,j,l)-LOx_HO2_O_c(t,z,j,l)
      Sum_D_HO2_O3(j)=Sum_D_HO2_O3(j)+LOx_HO2_O3_x(t,z,j,l)-LOx_HO2_O3_c(t,z,j,l)
      Sum_D_H_O3(j)=Sum_D_H_O3(j)+LOx_H_O3_x(t,z,j,l)-LOx_H_O3_c(t,z,j,l)
      Sum_D_OH_O(j)=Sum_D_OH_O(j)+LOx_OH_O_x(t,z,j,l)-LOx_OH_O_c(t,z,j,l)
      Sum_D_O_O(j)=Sum_D_O_O(j)+LOx_O_O_x(t,z,j,l)-LOx_O_O_c(t,z,j,l)
    end do
  end do
end do

do j=0,dims(2)-1
  Sum_loss=Sum_D_O3_O(j)+ Sum_D_NO2_O(j)+ Sum_ClO_O(j)+ Sum_HO2_O(j)+ Sum_HO2_O3(j)+ Sum_H_O3(j)+ Sum_OH_O(j)+ Sum_O_O(j)

  Area=cos(dble2flt(lat(j))*pi/180.)

  Sum_O3_O=Sum_O3_O+Sum_D_O3_O(j)*Sum_D_O3(j)*Area/Sum_loss

  Sum_NO2_O=Sum_NO2_O+Sum_D_NO2_O(j)*Sum_D_O3(j)*Area/Sum_loss

  Sum_ClO_O=Sum_ClO_O+Sum_D_ClO_O(j)*Sum_D_O3(j)*Area/Sum_loss

  Sum_HO2_O=Sum_HO2_O+Sum_D_HO2_O(j)*Sum_D_O3(j)*Area/Sum_loss

  Sum_HO2_O3=Sum_HO2_O3+Sum_D_HO2_O3(j)*Sum_D_O3(j)*Area/Sum_loss

  Sum_H_O3=Sum_H_O3+Sum_D_H_O3(j)*Sum_D_O3(j)*Area/Sum_loss

  Sum_OH_O=Sum_OH_O+Sum_D_OH_O(j)*Sum_D_O3(j)*Area/Sum_loss

  Sum_O_O=Sum_O_O+Sum_D_O_O(j)*Sum_D_O3(j)*Area/Sum_loss
  
  Sum_A=Sum_A+Area

  Sum_D_O3_glob=Sum_D_O3_glob+Sum_D_O3(j)*Area
  
end do
Sum_O3_O=Sum_O3_O/Sum_A
Sum_NO2_O=Sum_NO2_O/Sum_A
Sum_ClO_O=Sum_ClO_O/Sum_A
Sum_HO2_O=Sum_HO2_O/Sum_A
Sum_HO2_O3=Sum_HO2_O3/Sum_A
Sum_H_O3=Sum_H_O3/Sum_A
Sum_OH_O=Sum_OH_O/Sum_A
Sum_O_O=Sum_O_O/Sum_A
Sum_D_O3_glob=Sum_D_O3_glob/Sum_A

print(Sum_A)
print(Sum_O3_O)
print(Sum_NO2_O)
print(Sum_ClO_O)
print(Sum_HO2_O)
print(Sum_HO2_O3)
print(Sum_H_O3)
print(Sum_OH_O)
print(Sum_O_O)
print(Sum_D_O3_glob)
