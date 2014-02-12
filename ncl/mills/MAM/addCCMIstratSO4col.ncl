load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

  ccmifile="/glade/p/work/mmills/inputdata/CESM_1991_1995_sad_V2_sum_gw_SO2_TROPZ.nc"
  a=addfile(ccmifile,"w")

  tropz=a->TROP_Z
  so4=a->H2SO4_mass ; molecules/cm3
  factor=98.078479e3/6.022e23 ; molecules/cm3 -> kg/m3
  so4=so4*factor ; kg/m3
  nTimes=dimsizes(so4&time)
  nLats=dimsizes(so4&lat)
  
  so4_strat=a->H2SO4_col ; get metadata and array structure
  so4_strat@long_name="Stratospheric sulfate mass"
  so4_strat@units="kg/m2"
  
  do i = 0,nTimes-1
    print("month "+i)
    do l = 0,nLats-1
      tz=tropz(i,l)
      so4_strat(i,l)=dim_sum(so4(i, {tz:}, l))
    end do
  end do
  so4_strat(:,:)=so4_strat(:,:)*500.e0 ; 0.5 km altitude intervals
  
  a->so4_strat=so4_strat
  
  delete(a)
  
