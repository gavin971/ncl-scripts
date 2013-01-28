load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

root="/data/ccsm4/"

caseid="b40.1850.2deg.wset.002"
a=addfile(root+caseid+"/proc/"+caseid+".surfvars.nc","r")
wsetPI=a->surfvars
delete(a)

caseid="b40.1850.track1.1deg.006"
a=addfile(root+caseid+"/proc/"+caseid+".surfvars.nc","r")
camPI=a->surfvars
delete(a)

caseid="b40.1955-2005.2deg.wset.avg456"
a=addfile(root+caseid+"/proc/"+caseid+".surfvars.nc","r")
wset20th=a->surfvars
delete(a)

caseid="b40.20th.track1.1deg.avg567"
a=addfile(root+caseid+"/proc/"+caseid+".surfvars.nc","r")
cam20th=a->surfvars
delete(a)


root="/data/waccm/"

caseid="b40.1850.track1.2deg.wcm.007"
a=addfile(root+caseid+"/proc/"+caseid+".surfvars.nc","r")
wcmPI=a->surfvars
delete(a)

caseid="b40.1955-2005.2deg.wcm.avg"
a=addfile(root+caseid+"/proc/"+caseid+".surfvars.nc","r")
wcm20th=a->surfvars
latbounds=a->latbounds
delete(a)

do i = 0,3 ; loop over latitude bounds
  print(" ")
  print("Latitude bounds: "+latbounds(i,0)+" to "+latbounds(i,1))
  
  table=(/" WACCM& "," CCSM4& "," CCSM4-WSET& "/)
  do j=0,4
    diff1=wcm20th(i,j,0)-wcmPI(i,j,0)
    if (abs(diff1).gt.abs(2*wcmPI(i,j,1))) then
      table(0)=table(0)+"\textbf{"+sprintf("%5.2f", diff1)+"}& "
    else
      table(0)=table(0)+sprintf("%5.2f", diff1)+"& "
    end if
    
    diff1=cam20th(i,j,0)-camPI(i,j,0)
    if (abs(diff1).gt.abs(2*camPI(i,j,1))) then
      table(1)=table(1)+"\textbf{"+sprintf("%5.2f", diff1)+"}& "
    else
      table(1)=table(1)+sprintf("%5.2f", diff1)+"& "
    end if
    
    diff1=wset20th(i,j,0)-wsetPI(i,j,0)
    if (abs(diff1).gt.abs(2*wsetPI(i,j,1))) then
      table(2)=table(2)+"\textbf{"+sprintf("%5.2f", diff1)+"}& "
    else
      table(2)=table(2)+sprintf("%5.2f", diff1)+"& "
    end if
  end do
  print((/table/))
  delete(table)
  
end do
