function dnum = parse_time_units(time,units)
% datenum = parse_time_units(time,units)
%
% Parse a CF-compliant UNITS string from netcdf and apply the necessary scale 
% factor and offset to convert TIME to a datenum
%
% John Wilkin - July 2020

switch units(1:5)
  case 'milli'
    fac = 86400e3;
  case 'secon'
    fac = 86400;
  case 'minut'
    fac = 3600;
  case 'hours'
    fac = 24;
  case 'days '
    fac = 1;
end
% convert to datenum
tsince = units(6+strfind(units,'since'):end);
tsince = strrep(tsince,'UTC','');
dnum = time/fac + datenum(tsince);
