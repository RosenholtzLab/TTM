function strc = loadStructureFile(fname)
% strc = loadStructureFile(filename)
%
% Loads structure information from filename, used to deal with parameter
% and probably statistical structures saved during texture synthesis
% 
% filename: a file in the form
%     fieldname fieldvalue
%     where 
%       fieldname  is any string
%       fieldvalue is any number or string or matlab expression
%     strc.fieldname will have the given fieldvalue 
%
% Copyright (C) 2019  Ruth Rosenholtz, Alvin Raj, & Lavanya Sharan
% See COPYING.txt for license details

fid = fopen(fname);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    myline = strtrim(tline);
    if (numel(myline) == 0 || strcmp(myline(1),'%') == 1) 
        % is a comment or empty line, ignore
    else
        [first, rest] = strtok(tline);        
        rest = strtrim(rest); %remove leading and trailing whitespace
        try
            eval(sprintf('%s = [ %s ];', first, rest));
            eval(sprintf('strc.%s = [ %s ];', first, rest));
        catch
            % in case of problems, treat it as a string?
            eval(sprintf('%s = [ ''%s'' ];', first, rest)); 
            eval(sprintf('strc.%s = [ ''%s'' ];', first, rest));
        end
        
    end
end
fclose(fid);