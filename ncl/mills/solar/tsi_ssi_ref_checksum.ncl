ncf = addfile("spectra_1610_2009a_31Mar10.nc","r")

ssi_ref=ncf->ssi_ref
tsi_ref=ncf->tsi_ref

band_width=ncf->band_width
ssi_sum = ssi_ref # band_width
ssi_sum = ssi_sum/1000.d0 ; Convert mW to W

print(tsi_ref)
print(ssi_sum)
