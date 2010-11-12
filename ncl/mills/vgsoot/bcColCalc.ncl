
basedir="/Volumes/Data/Models/waccm/run/carma-soot/"

caseid1="bcrckt03"
rootdir1=basedir+caseid1+"/h0zm.annual/"
file1=caseid1+".cam2.h0zm.aer.annual.nc"

caseid2="bccntl02"
rootdir2=basedir+caseid2+"/h0zm.annual/"
file2=caseid2+".cam2.h0zm.aer.annual.nc"

bcrckt=addfile(rootdir1+file1,"r")
bccntl=addfile(rootdir2+file2,"r")

  lev=bcrckt->lev
  gw=bcrckt->gw
  
  MBCPHO_V=bcrckt->MBCPHO_V ; units=Kg/m^2
  MBCPHI_V=bcrckt->MBCPHI_V
  MOCPHO_V=bcrckt->MOCPHO_V
  MOCPHI_V=bcrckt->MOCPHI_V
  MBC_V_rckt=MBCPHO_V ; get attributes
  MBC_V_rckt=MBC_V_rckt+MBCPHI_V+MOCPHO_V+MOCPHI_V
  
  MBCPHO_V=bccntl->MBCPHO_V
  MBCPHI_V=bccntl->MBCPHI_V
  MOCPHO_V=bccntl->MOCPHO_V
  MOCPHI_V=bccntl->MOCPHI_V
  MBC_V_cntl=MBCPHO_V ; get attributes
  MBC_V_cntl=MBC_V_cntl+MBCPHI_V+MOCPHO_V+MOCPHI_V
  
  bounds=(/(/-91,0/),(/0,91/),(/-91,91/)/)
  titles=(/"Southern Hemisphere","Northern Hemisphere","Whole Earth"/)
  do i = 0,2

    masstop_rckt=MBC_V_rckt(0,{(/bounds(i,0)/):(/bounds(i,1)/)})#gw({(/bounds(i,0)/):(/bounds(i,1)/)})
    masstop_cntl=MBC_V_cntl(0,{(/bounds(i,0)/):(/bounds(i,1)/)})#gw({(/bounds(i,0)/):(/bounds(i,1)/)})
    masstrop_rckt=MBC_V_rckt({100},{(/bounds(i,0)/):(/bounds(i,1)/)})#gw({(/bounds(i,0)/):(/bounds(i,1)/)})
    masstrop_cntl=MBC_V_cntl({100},{(/bounds(i,0)/):(/bounds(i,1)/)})#gw({(/bounds(i,0)/):(/bounds(i,1)/)})

    masstop_rckt=masstop_rckt*510065600*1e6 ; Surface area of the Earth (m2)
    masstop_cntl=masstop_cntl*510065600*1e6 ; Surface area of the Earth (m2)
    masstrop_rckt=masstrop_rckt*510065600*1e6 ; Surface area of the Earth (m2)
    masstrop_cntl=masstrop_cntl*510065600*1e6 ; Surface area of the Earth (m2)

    massstrat_rckt=masstop_rckt-masstrop_rckt
    massstrat_cntl=masstop_cntl-masstrop_cntl
    massstrat_enhance=(massstrat_rckt/massstrat_cntl-1.0)*100.0
    massstrat_diff=massstrat_rckt-massstrat_cntl

    print(" ")  
    print("====================================================") 
    print((/titles(i)/))
    print("====================================================")  
    print("Rocket total BC:                      "+masstop_rckt+" kg")
    print("Rocket tropospheric BC   (p>100 hPa): "+masstrop_rckt+" kg")
    print("Rocket stratospheric BC  (p<100 hPa): "+massstrat_rckt+" kg")
    print("====================================================")  
    print("Control total BC:                     "+masstop_cntl+" kg")
    print("Control tropospheric BC  (p>100 hPa): "+masstrop_cntl+" kg")
    print("Control stratospheric BC (p<100 hPa): "+massstrat_cntl+" kg")
    print("====================================================") 
    print("Rocket stratospheric BC enhancement:  "+massstrat_enhance+"%")
    print("Rocket-Control stratospheric BC:      "+massstrat_diff+" kg")
    print("====================================================")  

  end do
