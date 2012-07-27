datapath="/Volumes/MillsPassport/cesm/case/"

ccsmPI="b40.1850.track1.1deg.006"
ccsm20th="b40.20th.track1.1deg.avg5678912"

waccmPI="b40.1850.track1.2deg.wcm.007"
waccm20th="b40.1955-2005.2deg.wcm.avg234"
waccm20th_cases="b40.1955-2005.2deg.wcm.00"+ispan(2,4,1)

subdir="/proc/tseries/monthly/"

;print("opening "+ccsmPI+".cam2.h0gm.TS.080001-099912.nc")
ccsmPI_TSfile=addfile(datapath+ccsmPI+subdir+ccsmPI+".cam2.h0gm.TS.080001-099912.nc","r")
;print("opening "+ccsmPI+".cam2.h0nh.TS.080001-099912.nc")
ccsmPI_TSnhfile=addfile(datapath+ccsmPI+subdir+ccsmPI+".cam2.h0nh.TS.080001-099912.nc","r")
;print("opening "+ccsmPI+".cam2.h0sh.TS.080001-099912.nc")
ccsmPI_TSshfile=addfile(datapath+ccsmPI+subdir+ccsmPI+".cam2.h0sh.TS.080001-099912.nc","r")
;print("opening "+ccsmPI+".cam2.h0gm.FSNT.080001-099912.nc")
ccsmPI_FSfile=addfile(datapath+ccsmPI+subdir+ccsmPI+".cam2.h0gm.FSNT.080001-099912.nc","r")
;print("opening "+ccsmPI+".cam2.h0gm.FLNT.080001-099912.nc")
ccsmPI_FLfile=addfile(datapath+ccsmPI+subdir+ccsmPI+".cam2.h0gm.FLNT.080001-099912.nc","r")

ccsmPI_TS=ccsmPI_TSfile->TS
ccsmPI_TSnh=ccsmPI_TSnhfile->TS
ccsmPI_TSsh=ccsmPI_TSshfile->TS
ccsmPI_date=ccsmPI_TSfile->date
ccsmPI_FS=ccsmPI_FSfile->FSNT
ccsmPI_FL=ccsmPI_FLfile->FLNT
ccsmPI_RESTOM=ccsmPI_FS
ccsmPI_RESTOM=(/ccsmPI_FS-ccsmPI_FL/)
print("ccsm PI 080001-099912: TSga, TSnh, TSsh, RESTOM, FSNT, FLNT")
print(avg(ccsmPI_TS)+" "+avg(ccsmPI_TSnh)+" "+avg(ccsmPI_TSsh)+" "+avg(ccsmPI_RESTOM)+" "+avg(ccsmPI_FS)+" "+avg(ccsmPI_FL))

;print("opening "+ccsm20th+".cam2.h0gm.TS.185001-200512.nc")
ccsm20th_TSfile=addfile(datapath+ccsm20th+subdir+ccsm20th+".cam2.h0gm.TS.185001-200512.nc","r")
;print("opening "+ccsm20th+".cam2.h0nh.TS.185001-200512.nc")
ccsm20th_TSnhfile=addfile(datapath+ccsm20th+subdir+ccsm20th+".cam2.h0nh.TS.185001-200512.nc","r")
;print("opening "+ccsm20th+".cam2.h0sh.TS.185001-200512.nc")
ccsm20th_TSshfile=addfile(datapath+ccsm20th+subdir+ccsm20th+".cam2.h0sh.TS.185001-200512.nc","r")
;print("opening "+ccsm20th+".cam2.h0gm.FSNT.185001-200512.nc")
ccsm20th_FSfile=addfile(datapath+ccsm20th+subdir+ccsm20th+".cam2.h0gm.FSNT.185001-200512.nc","r")
;print("opening "+ccsm20th+".cam2.h0gm.FNLT.185001-200512.nc")
ccsm20th_FLfile=addfile(datapath+ccsm20th+subdir+ccsm20th+".cam2.h0gm.FLNT.185001-200512.nc","r")

ccsm20th_TS=ccsm20th_TSfile->TS
ccsm20th_TSnh=ccsm20th_TSnhfile->TS
ccsm20th_TSsh=ccsm20th_TSshfile->TS
ccsm20th_date=ccsm20th_TSfile->date
ccsm20th_FS=ccsm20th_FSfile->FSNT
ccsm20th_FL=ccsm20th_FLfile->FLNT
ccsm20th_RESTOM=ccsm20th_FS
ccsm20th_RESTOM=(/ccsm20th_FS-ccsm20th_FL/)
indrange=ind(ccsm20th_date.gt.19860101.and.ccsm20th_date.le.20060101)
print(" ")
print("ccsm 1986-2005: TSga, TSnh, TSsh, RESTOM, FSNT, FLNT")
print(avg(ccsm20th_TS(indrange))+" "+avg(ccsm20th_TSnh(indrange))+" "+avg(ccsm20th_TSsh(indrange))+" "\
+avg(ccsm20th_RESTOM(indrange))+" "+avg(ccsm20th_FS(indrange))+" "+avg(ccsm20th_FL(indrange)))
print("TS diff "+(/avg(ccsm20th_TS(indrange))-avg(ccsmPI_TS)/))
print("FS diff "+(/avg(ccsm20th_FS(indrange))-avg(ccsmPI_FS)/))
print("FL diff "+(/avg(ccsm20th_FL(indrange))-avg(ccsmPI_FL)/))
delete(indrange)

;print("opening "+waccmPI+".cam2.h0gm.TS.009601-029512.nc")
waccmPI_TSfile=addfile(datapath+waccmPI+subdir+waccmPI+".cam2.h0gm.TS.009601-029512.nc","r")
;print("opening "+waccmPI+".cam2.h0nh.TS.009601-029512.nc")
waccmPI_TSnhfile=addfile(datapath+waccmPI+subdir+waccmPI+".cam2.h0nh.TS.009601-029512.nc","r")
;print("opening "+waccmPI+".cam2.h0sh.TS.009601-029512.nc")
waccmPI_TSshfile=addfile(datapath+waccmPI+subdir+waccmPI+".cam2.h0sh.TS.009601-029512.nc","r")
;print("opening "+waccmPI+".cam2.h0gm.FSNT.009601-029512.nc")
waccmPI_FSfile=addfile(datapath+waccmPI+subdir+waccmPI+".cam2.h0gm.FSNT.009601-029512.nc","r")
;print("opening "+waccmPI+".cam2.h0gm.FLNT.009601-029512.nc")
waccmPI_FLfile=addfile(datapath+waccmPI+subdir+waccmPI+".cam2.h0gm.FLNT.009601-029512.nc","r")

waccmPI_TS=waccmPI_TSfile->TS
waccmPI_TSnh=waccmPI_TSnhfile->TS
waccmPI_TSsh=waccmPI_TSshfile->TS
waccmPI_date=waccmPI_TSfile->date
waccmPI_FS=waccmPI_FSfile->FSNT
waccmPI_FL=waccmPI_FLfile->FLNT
waccmPI_RESTOM=waccmPI_FS
waccmPI_RESTOM=(/waccmPI_FS-waccmPI_FL/)
print(" ")
print("waccm PI 009601-029512: TSga, TSnh, TSsh, RESTOM, FSNT, FLNT")
print(avg(waccmPI_TS)+" "+avg(waccmPI_TSnh)+" "+avg(waccmPI_TSsh)+" "+avg(waccmPI_RESTOM)+" "+\
 avg(waccmPI_FS)+" "+avg(waccmPI_FL))

;print("opening "+waccm20th+".cam2.h0gm.TS.195501-200512.nc")
waccm20th_TSfile=addfile(datapath+waccm20th+subdir+waccm20th+".cam2.h0gm.TS.195501-200512.nc","r")
;print("opening "+waccm20th+".cam2.h0nh.TS.195501-200512.nc")
waccm20th_TSnhfile=addfile(datapath+waccm20th+subdir+waccm20th+".cam2.h0nh.TS.195501-200512.nc","r")
;print("opening "+waccm20th+".cam2.h0sh.TS.195501-200512.nc")
waccm20th_TSshfile=addfile(datapath+waccm20th+subdir+waccm20th+".cam2.h0sh.TS.195501-200512.nc","r")
;print("opening "+waccm20th+".cam2.h0gm.FSNT.195501-200512.nc")
waccm20th_FSfile=addfile(datapath+waccm20th+subdir+waccm20th+".cam2.h0gm.FSNT.195501-200512.nc","r")
;print("opening "+waccm20th+".cam2.h0gm.FLNT.195501-200512.nc")
waccm20th_FLfile=addfile(datapath+waccm20th+subdir+waccm20th+".cam2.h0gm.FLNT.195501-200512.nc","r")

waccm20th_TS=waccm20th_TSfile->TS
waccm20th_TSnh=waccm20th_TSnhfile->TS
waccm20th_TSsh=waccm20th_TSshfile->TS
waccm20th_date=waccm20th_TSfile->date
waccm20th_FS=waccm20th_FSfile->FSNT
waccm20th_FL=waccm20th_FLfile->FLNT
waccm20th_RESTOM=waccm20th_FS
waccm20th_RESTOM=(/waccm20th_FS-waccm20th_FL/)

indrange=ind(waccm20th_date.gt.19860101.and.waccm20th_date.le.20060101)
print(" ")
print("WACCM 1986-2005: TSga, TSnh, TSsh, RESTOM, FSNT, FLNT")
print(avg(waccm20th_TS(indrange))+" "+avg(waccm20th_TSnh(indrange))+" "+avg(waccm20th_TSsh(indrange))+\
" "+avg(waccm20th_RESTOM(indrange))+" "+avg(waccm20th_FS(indrange))+" "+avg(waccm20th_FL(indrange)))
print("TS diff "+(/avg(waccm20th_TS(indrange))-avg(waccmPI_TS)/))
print("FS diff "+(/avg(waccm20th_FS(indrange))-avg(waccmPI_FS)/))
print("FL diff "+(/avg(waccm20th_FL(indrange))-avg(waccmPI_FL)/))

do i =0,2
  delete(waccm20th_TSfile)
  delete(waccm20th_TS)
  waccm20th_TSfile=addfile(datapath+waccm20th_cases(i)+subdir+waccm20th_cases(i)+".cam2.h0gm.TS.195501-200512.nc","r")
  waccm20th_TS=waccm20th_TSfile->TS
  print(waccm20th_cases(i)+": TS diff "+(/avg(waccm20th_TS(indrange))-avg(waccmPI_TS)/))
end do
