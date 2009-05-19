function getDataPath ()
local dataPath
begin

  if (isfilepresent("/Volumes/ptmp")) then
    dataPath="/Volumes/ptmp/mills/Data/"
  else
    if (isfilepresent("/ptmp")) then
      dataPath="/ptmp/mills/Data/"
    else
      if (isfilepresent("~/Data")) then
        dataPath="~/Data/"
      end if
    end if
  end if

  return(dataPath)
end
