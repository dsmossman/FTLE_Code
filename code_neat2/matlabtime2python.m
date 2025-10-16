function [pythontime,date_time]=matlabtime2python(dnum)

date_time = datetime(dnum,'ConvertFrom','datenum');
pythontime = seconds(date_time - datetime(1970,1,1));

return