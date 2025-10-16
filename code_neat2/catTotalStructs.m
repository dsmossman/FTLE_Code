function [TUVcat,goodCount] = catTotalStructs(f,v,ef,of)

% f = cell array containing names of mat files with TUV structs.
% v = name of the struct in the mat file to concat
% ef = concat error estimates - defaults to false
% of = concat other vars - defaults to true

% TUVcat = output concat'ed tuv struct
% cnt = number of successfully loaded and concat'ed structs.

if ~exist('ef','var') || isempty(ef)
    ef = false;
end
if ~exist('of','var') || isempty(ef)
    of = true;
end
if ~exist('v','var')  || isempty(v)
    v = 'TUV';
end

numFiles = length(f);
goodCount = 0;
for i = 1:numFiles
    try
        % Load total hour:  assume it has the struct TUV as the final total current
        % data product.  totals or oma should have TUV if using
        % HFRPdriver_Totals_OMA.m
        % Assuming if I can load it it's good.  What if the file exists but
        % is empty?  I think it will get counted as a good file.
        data=load(f{i},v);
        fprintf('Loading file %u of %u\n',i,numFiles);
   catch
        fprintf('Can''t load %s ... skipping\n',f{i});
        continue;  % Skip rest of for loop
    end
        
    TUV = data.(v);
    
    TUV.OtherTemporalVars.OMA_NumRadialsPerSite = [];
    
    goodCount = goodCount + 1;
    
    if goodCount == 1
        TUVcat = TUV;
    end
    % Concat temporally - do a try catch because some temporalConcatTUV's
    % don't work - the problem is with the
    % OtherTemporalVars.OMA_NumRadialsPerSite field.  If the number of
    % sites changes over time, this field is a vector whos number of
    % elements changes and temporalConcatTUV will error.  So get rid of it
    % because it appears to be informational, not necessary.  But can't set
    % the catOthers to false because the particle analysis needs the
    % OtherTemporalVars.OMA_alpha's.
    % This was set to empty above, so shouldn't hit the catch part now. 
    if goodCount > 1
        try
            TUVcat = temporalConcatTUV(TUVcat,TUV,ef,of);
        catch
            fprintf('Can''t concat %s ... skipping\n',f{i});
            goodCount = goodCount - 1;
        end
    end

    clear data
end
