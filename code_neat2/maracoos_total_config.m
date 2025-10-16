function conf = maracoos_total_config(domain)

compType = computer;

if ~isempty(strmatch('MACI64', compType))
     root = '/Volumes';
else
     root = '/home';
end

%% radial configs
conf.Radials.BaseDir= [root '/codaradm/data/radials/'];
conf.Radials.MonthFlag=false;
conf.Radials.TypeFlag=false;
conf.Radials.MaskDir=[root '/codaradm/data/mask_files/MARACOOS_6kmMask.txt'];

if exist('domain')
    if strcmp(domain,'MARACOOS')
        conf.Radials.Sites= {'NAUS','NAUS','NANT','NANT','MVCO','MVCO','BLCK','BLCK'};
        conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='MARACOOS';
        
        conf.Totals.Type='TUV';
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:5:30;
        conf.HourPlot.axisLims=[-77 -67 33 44];
        conf.AnimPlot.axisLims=[-77 -67 33 44];
        conf.HourPlot.DomainName='MARA';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/codaradm/data/coast_files/MARA_coast.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        conf.HourPlot.Print=false;
       
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/maracoos/oi/mat/5MHz/'];
        conf.OI.FilePrefix = 'tuv_oi_MARA_';
       % removed this path on 5/20/22
       % conf.OI.BaseDir=[ root '/codaradm/data/totals/realtime/maracoos/oi/mat/5MHz/'];
        %conf.OI.FilePrefix='hfr_rtv_midatl_6km_oi_maracoos_';


        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        
        conf.gridspacing=3;
        
        addpath /home/hroarty/codar/MARACOOS/Animation_Total/
        
%         wp1=[40+30/60 -71-20/60];% removed 20161004 HJR   
        wp1=[40+30/60 -70-00/60];
        resolution=6;
        range=60;
        [conf.Animation.wp]=release_point_generation_matrix(wp1,resolution,range);
        conf.Animation.PrintPath='/www/home/hroarty/public_html/maracoos/animations/MARACOOS_North/';
        conf.Animation.buffer=15/60;
        
    elseif strcmp(domain,'MARACOOS_Central')
        conf.Radials.Sites= {'NAUS','NAUS','NANT','NANT','MVCO','MVCO','BLCK','BLCK'};
        conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='MARACOOS_Central';
        
        conf.Totals.Type='TUV';
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:5:30;
        conf.HourPlot.axisLims=[-75 -73 38.25 40.75];
        conf.AnimPlot.axisLims=[-75 -73 38.25 40.75];
        conf.HourPlot.DomainName='MARA';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/codaradm/data/coast_files/MARA_coast.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        conf.HourPlot.Print=false;
       
        
        % conf.OI.BaseDir=[root '/codaradm/data/totals/maracoos/oi/mat/5MHz/'];
        %  conf.OI.FilePrefix = 'tuv_oi_MARA_';
        % removed the above path on 5/20/22
         conf.OI.BaseDir=[ root '/codaradm/data/totals/realtime/maracoos/oi/mat/5MHz/'];
         conf.OI.FilePrefix='hfr_rtv_midatl_6km_oi_maracoos_';
        
        
        
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        
        conf.gridspacing=1;
        
        addpath /home/hroarty/codar/MARACOOS/Animation_Total/
        
%         wp1=[40+30/60 -71-20/60];% removed 20161004 HJR   
        wp1=[40+30/60 -70-00/60];
        resolution=6;
        range=60;
        [conf.Animation.wp]=release_point_generation_matrix(wp1,resolution,range);
        conf.Animation.PrintPath='/www/home/hroarty/public_html/maracoos/animations/MARACOOS_Central/';
        conf.Animation.buffer=15/60;
        
    elseif strcmp(domain,'Thirteen')
        conf.Radials.Sites= {'BELM','BELM','BRNT','BRNT','BRMR','BRMR','RATH','RATH','WOOD','WOOD'};
        conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='Thirteen';
        
        conf.Totals.Type='TUV';
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:5:30;
        conf.HourPlot.axisLims=[-75 -73 38.25 40.75];
        conf.HourPlot.DomainName='BPU';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/hroarty/data/coast_files/13MHz_NJ.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/maracoos/oi/mat/13MHz/'];
        conf.OI.FilePrefix=['tuv_oi_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        
         conf.gridspacing=2;
         
        addpath /home/hroarty/codar/MARACOOS/Animation_Total/
        wp1=[39+30/60 -74-00/60];   
        resolution=4;
        range=40;
        [conf.Animation.wp]=release_point_generation_matrix(wp1,resolution,range);
        conf.Animation.PrintPath='/www/home/hroarty/public_html/maracoos/animations/Thirteen/';
        conf.Animation.buffer=10/60;
        
     elseif strcmp(domain,'ThirteenQC')
        conf.Radials.Sites= {'BELM','BELM','BRNT','BRNT','BRMR','BRMR','RATH','RATH','WOOD','WOOD'};
        conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='ThirteenQC';
        
        conf.Totals.Type='TUVosn';
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:5:30;
        conf.HourPlot.axisLims=[-75 -73 38.25 40.75];
        conf.HourPlot.DomainName='BPU';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/hroarty/data/coast_files/13MHz_NJ.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/hroarty/data/realtime/totals/maracoos/oi/mat/13MHz/'];
        conf.OI.FilePrefix=['tuv_oi_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        
         conf.gridspacing=2;
        
    elseif strcmp(domain,'NYHarbor')
        conf.Radials.Sites= {'PORT','PORT','SILD','SILD'};
        conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm'};
        conf.Region='NYHarbor';
        
        conf.Totals.Type='TUV';
        
        conf.Radials.MaskDir=[root '/codaradm/data/mask_files/25MHz_1kmMask.txt'];
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:5:30;
        conf.HourPlot.axisLims=[-74.33 -73.75  40.33 40.66];
        conf.AnimPlot.axisLims=[-74.33 -73.75  40.33 40.66];
        conf.HourPlot.DomainName='MARASR';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/codaradm/data/coast_files/NYH_coast.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/maracoos/oi/mat/25MHz/'];
        conf.OI.FilePrefix=['tuv_oi_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        conf.OI.Mask=[root '/hroarty/data/mask_files/25MHz_Mask.txt'];
        
        conf.gridspacing=1;
       
        addpath /home/hroarty/codar/MARACOOS/Animation_Total/
        [conf.Animation.wp]=release_point_generation_matrix_NYH;
        conf.Animation.PrintPath='/www/home/hroarty/public_html/maracoos/animations/NYH/';
        conf.Animation.buffer=5/60;
    
    elseif strcmp(domain,'PuertoRico')
        conf.Radials.Sites= {'CDDO','CDDO','FURA','FURA'};
        conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm'};
        conf.Region='PuertoRico';
        
        conf.Totals.Type='TUV';
        
        conf.Radials.MaskDir=[root '/codaradm/data/mask_files/PR_coast.mask'];
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:10:50;
        conf.HourPlot.axisLims=[-68.30 -67 17.5 18.75];
        conf.AnimPlot.axisLims=[-68.30 -67 17.5 18.75];
        conf.HourPlot.DomainName='CARA';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/hroarty/data/coast_files/UPRM_coast.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/caracoos/oi/mat/13MHz/measured/'];
        conf.OI.FilePrefix=['tuv_OI_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        conf.OI.Mask=[root '/codaradm/data/mask_files/PR_coast.mask'];
        
        conf.gridspacing=2;
        
        addpath /home/hroarty/codar/MARACOOS/Animation_Total/
        wp1=[17+45/60 -68-00/60];
        resolution=8;
        range=80;
        [conf.Animation.wp]=release_point_generation_matrix(wp1,resolution,range);
        conf.Animation.PrintPath='/www/home/hroarty/public_html/maracoos/animations/PR/';
        conf.Animation.buffer=15/60;
        
        
    elseif strcmp(domain,'PLDP')
        conf.Radials.Sites= {'JOUB','JOUB','WAUW','WAUW','PALM','PALM'};
        conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='PLDP';
        
        conf.Totals.Type='TUV';
        
        conf.Radials.MaskDir=[root '/hroarty/data/mask_files/Antarctica_Coast_Mask.txt'];
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:10:50;
        conf.HourPlot.axisLims=[-65.0 -63.5 -65.25 -64.5 ];% lon lat
        conf.AnimPlot.axisLims=[-64.5 -63.833 -65.1 -64.75];
        conf.HourPlot.DomainName='PLDP';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/hroarty/data/coast_files/Antarctica4.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/pldp/25MHz/0.5km/oi/mat'];
        conf.OI.FilePrefix=['tuv_oi_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        conf.OI.Mask=[root '/codaradm/data/mask_files/PLDP_coast.mask'];
        
        conf.gridspacing=2;
        
        addpath /home/hroarty/codar/MARACOOS/Animation_Total/
        [conf.Animation.wp]=release_point_generation_matrix_PLDP;
        conf.Animation.PrintPath='/www/home/hroarty/public_html/maracoos/animations/PLDP/';
        conf.Animation.buffer=15/60;
        
    elseif strcmp(domain,'DelawareBay')
%         conf.Radials.Sites= {'JOUB','JOUB','WAUW','WAUW','PALM','PALM'};
%         conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='DelawareBay';
        
        conf.Totals.Type='TUV';
        
        conf.Radials.MaskDir=[root '/hroarty/data/mask_files/Antarctica_Coast_Mask.txt'];
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:10:50;
        conf.HourPlot.axisLims=[-75-20/60 -74-40/60 38+30/60 39+12/60 ];% lon lat
        conf.AnimPlot.axisLims=[-75-20/60 -74-40/60 38+30/60 39+12/60];
        conf.HourPlot.DomainName='MARASR';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/codaradm/data/coast_files/MARA_coast.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/maracoos/oi/mat/25MHz/'];
        conf.OI.FilePrefix=['tuv_oi_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        conf.OI.Mask=[root '/hroarty/data/mask_files/25MHz_Mask_Delaware_Bay.txt'];
        
        conf.gridspacing=1;

        addpath /home/hroarty/codar/MARACOOS/Animation_Total/
        wp1=[38+40/60 -75-12/60];
        resolution=3;
        range=40;
        [conf.Animation.wp]=release_point_generation_matrix(wp1,resolution,range);
        conf.Animation.PrintPath='/www/home/hroarty/public_html/maracoos/animations/DelawareBay/';
        conf.Animation.buffer=10/60;
        
    elseif strcmp(domain,'BlockIslandSound')
%         conf.Radials.Sites= {'JOUB','JOUB','WAUW','WAUW','PALM','PALM'};
%         conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='BlockIslandSound';
        
        conf.Totals.Type='TUV';
        
%         conf.Radials.MaskDir=[root '/hroarty/data/mask_files/Antarctica_Coast_Mask.txt'];
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:10:50;
        conf.HourPlot.axisLims=[-72-10/60 -71-10/60 40+40/60 41+30/60 ];% lon lat
        conf.AnimPlot.axisLims=[-72-10/60 -71-10/60 40+40/60 41+30/60];
        conf.HourPlot.DomainName='MARASR';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/codaradm/data/coast_files/MARA_coast.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/maracoos/oi/mat/25MHz/'];
        conf.OI.FilePrefix=['tuv_oi_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        conf.OI.Mask=[root '/hroarty/data/mask_files/25MHz_Mask.txt'];
        
        conf.gridspacing=2;
        
    elseif strcmp(domain,'LongIslandSound')
%         conf.Radials.Sites= {'JOUB','JOUB','WAUW','WAUW','PALM','PALM'};
%         conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='LongIslandSound';
        
        conf.Totals.Type='TUV';
        
%         conf.Radials.MaskDir=[root '/hroarty/data/mask_files/Antarctica_Coast_Mask.txt'];
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:10:50;
        conf.HourPlot.axisLims=[-73-50/60 -73-10/60 40+50/60 41+10/60 ];% lon lat
        conf.AnimPlot.axisLims=[-73-50/60 -73-10/60 40+50/60 41+10/60];
        conf.HourPlot.DomainName='MARASR';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/codaradm/data/coast_files/MARA_coast.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/maracoos/oi/mat/25MHz/'];
        conf.OI.FilePrefix=['tuv_oi_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        conf.OI.Mask=[root '/hroarty/data/mask_files/25MHz_Mask.txt'];
        
        conf.gridspacing=1;
        
    elseif strcmp(domain,'ChesBay')
%         conf.Radials.Sites= {'JOUB','JOUB','WAUW','WAUW','PALM','PALM'};
%         conf.Radials.Types= {'RDLi','RDLm','RDLi','RDLm','RDLi','RDLm'};
        conf.Region='ChesBay';
        
        conf.Totals.Type='TUV';
        
%         conf.Radials.MaskDir=[root '/hroarty/data/mask_files/Antarctica_Coast_Mask.txt'];
        
        conf.HourPlot.Type='OI';%% this can be OI or UWLS
        conf.HourPlot.VectorScale=0.004;
        conf.HourPlot.ColorMap='jet';
        conf.HourPlot.ColorTicks=0:10:50;
        conf.HourPlot.axisLims=[-76-30/60 -75-40/60 36+50/60 37+24/60 ];% lon lat
        conf.AnimPlot.axisLims=[-76-30/60 -75-40/60 36+50/60 37+24/60];
        conf.HourPlot.DomainName='MARASR';
        conf.HourPlot.Print=false;
        conf.HourPlot.coastFile=[root '/codaradm/data/coast_files/MARA_coast.mat'];
        conf.HourPlot.Projection='Mercator';
        conf.HourPlot.plotBasemap_xargs={'patch',[1 .69412 0.39216],'edgecolor','k'};
        
        
        conf.OI.BaseDir=[root '/codaradm/data/totals/maracoos/oi/mat/25MHz/'];
        conf.OI.FilePrefix=['tuv_oi_' conf.HourPlot.DomainName '_'];
        conf.OI.FileSuffix='.mat';
        conf.OI.MonthFlag=true;
        conf.OI.Mask=[root '/hroarty/data/mask_files/25MHz_Mask.txt'];
        
        conf.gridspacing=1;
        
    else
        error('Unknown Domain')
    end
end

conf.HourPlot.PrintPath= '/www/home/hroarty/public_html/maracoos/total_coverage/img/';
% conf.HourPlot.PrintPath= '/Volumes/www_hroarty/maracoos/total_coverage/img/';

conf.Totals.Mask2=conf.HourPlot.axisLims([1 3;2 3;2 4;1 4;1 3]);
conf.Totals.MaxTotSpeed=300;

conf.OI.cleanTotalsVarargin={{'OIuncert','Uerr',0.8},{'OIuncert','Verr',0.8}};

conf.CoveragePlot.ColorTicks=0:10:100;
conf.CoveragePlot.ColorMap='jet';






% conf.CoveragePlot.ColorTicks=0:10:100;
% conf.CoveragePlot.ColorMap='jet';

end



