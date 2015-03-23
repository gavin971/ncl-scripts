depjob=getenv("depjob")
caseid=getenv("caseid")
resubmitstr=getenv("resubmit")
widths=getenv("widths")

if ismissing(resubmitstr) then
  resubmit=0
else
  resubmit=stringtoint(resubmitstr)
end if


if ismissing(caseid) then
  print("FATAL: caseid environment variable not set")
  status_exit(1)
else
  if ismissing(depjob) then
    depstr=""
  else
    depstr=" -w "+depjob
  end if
  do i = 0, resubmit
    command="cd /glade/p/work/mmills/case/"+caseid+"/ ; bsub"+depstr+" < "+caseid+".run"
    print((/command/))
    a=systemfunc(command)
    print((/a/)) ;      Job <980757> is submitted to queue <regular>.
    depjob=str_get_field(a,2,"<>")
    print("depjob="+depjob)
    depstr=" -w "+depjob
    if .not.ismissing(widths) then
      print("widths="+widths)
      command="bsub -n 1 -q geyser -W 24:00 -P P93300007"+depstr+" ncl /glade/u/home/mmills/scripts/ncl/MAM/calcAerProp.ncl"
      print((/command/))
      a=systemfunc(command)
      print((/a/)) ;      Job <980757> is submitted to queue <regular>.
    end if
  end do
end if
