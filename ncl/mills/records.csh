#! /bin/csh

# Create annual, monthly and seasonal averages from an annual file.

# Annual
# d - n
ncra h0.$1.nc h0.$1.ave.nc

# Seasons
# djf, mam, jja, son
ncra -d time,1,90     h0.$1.nc h0.$1.1.nc
ncra -d time,91,181   h0.$1.nc h0.$1.2.nc
ncra -d time,182,273  h0.$1.nc h0.$1.3.nc
ncra -d time,274,364  h0.$1.nc h0.$1.4.nc
ncrcat h0.$1.?.nc h0.$1.sea.nc
rm h0.$1.?.nc

# Months
# d, j, f, m, a, m, j, j, a, s, o, n
ncra -d time,1,31     h0.$1.nc h0.$1.12.nc
ncra -d time,32,62    h0.$1.nc h0.$1.01.nc
ncra -d time,63,90    h0.$1.nc h0.$1.02.nc
ncra -d time,91,121   h0.$1.nc h0.$1.03.nc
ncra -d time,122,151  h0.$1.nc h0.$1.04.nc
ncra -d time,152,182  h0.$1.nc h0.$1.05.nc
ncra -d time,183,212  h0.$1.nc h0.$1.06.nc
ncra -d time,213,243  h0.$1.nc h0.$1.07.nc
ncra -d time,244,273  h0.$1.nc h0.$1.08.nc
ncra -d time,274,303  h0.$1.nc h0.$1.09.nc
ncra -d time,304,334  h0.$1.nc h0.$1.10.nc
ncra -d time,335,364  h0.$1.nc h0.$1.11.nc
ncrcat h0.$1.12.nc h0.$1.0?.nc h0.$1.10.nc h0.$1.11.nc h0.$1.mon.nc
rm h0.$1.??.nc
