; Program creates ssi_ref and tsi_ref values for average
; of years 1834-1867 (solar cycles 8-10)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

ncf = addfile("spectra_1610_2009a_31Mar10.nc","w")

ssi=ncf->ssi
tsi=ncf->tsi
date=ncf->date
time=ncf->time
print(date(224:257)+" "+tsi(224:257))

tsi_ref=dim_avg_Wrap(tsi(224:257))
print(tsi_ref)
;print()

ssi_ref=dim_avg_Wrap(ssi(wavelength|:,time|224:257))
print(ssi_ref(0:10))

time@bounds="time_bound"

time_bound=(/time(224),time(257)/)
time_bound@long_name = "time interval boundaries"
time_bound@units = time@units
time_bound!0="time_bound"

tsi_ref!0="average"
tsi_ref@long_name="average of tsi over time_bounds inclusive"
tsi_ref@time_op = "average"
;                tsi_ref:_FillValue = -9999. ;
;                tsi_ref:average_op_ncl = "dim_avg over dimension(s): time" ;

ssi_ref@long_name="average of ssi over time_bounds inclusive"
ssi_ref@time_op = "average"

ncf->time=time
ncf->time_bound=time_bound
ncf->tsi_ref=tsi_ref
ncf->ssi_ref=ssi_ref

