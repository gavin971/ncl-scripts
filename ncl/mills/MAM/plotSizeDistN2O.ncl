load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

dlog10Dp = 0.056268 ; bin width of Wilson FCAS measurements

; Bin geometric mean diameters for Wilson FCAS measurements:
DpMean = (/ 0.064015516,0.072870752,0.082950928,0.094425491,0.107487324, \
            0.122355994,0.139281441,0.158548176,0.18048007,0.205445792, \
            0.233865011,0.266215446,0.303040902,0.344960406,0.392678615, \
            0.446997661,0.508830634,0.579216934,0.659339737,0.750545889, \
            0.854368544,0.972552939,1.107085724,1.260228365,1.434555155, \
            1.632996487,1.858888114,2.116027222,2.40873626,2.741935601, \
            3.121226247 /)


sigmag=(/1.6, 1.6, 1.8/) ; geometric standard deviations for 3 MAM3 modes

rmode=
