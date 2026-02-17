if (isStrSubset(var,"T")) then
  opt         = True
  opt@bin_min = 254
  opt@bin_max = 308
end if
if (isStrSubset(var,"Q")) then
  opt         = True
  opt@bin_min = 0
  opt@bin_max = 18
end if

if (isStrSubset(var,"U")) then
  opt         = True
  opt@bin_min = -15
  opt@bin_max =  15
end if

if (isStrSubset(var,"V")) then
  opt         = True
  opt@bin_min = -15
  opt@bin_max =  15
end if

