function C = LoadXCATField(filename, SLICE, useInterpolated,single)
%% LoadXCAT_motion_fields

if nargin < 3
    useInterpolated=true;
end

if nargin < 1 || isempty(filename)
    % get MRXCAT output (.txt file)
    [filename, pathname] = uigetfile({'*.txt'},'','~/results/XCAT');
    filename = [pathname filename];
end

if nargin < 2 
    SLICE = 1;
end

basename = filename(1:end-5);

if single==0
    files = dir([basename '*.txt']); Nfiles = length(files);
    file_vec= 2:Nfiles+1;
else
    file_vec=single;
end


for n =  file_vec
    
    fprintf('\n n = %d ...', n);
    filename = [basename num2str(n) '.txt'];
    fid = fopen(filename,'rt');
    fgets(fid);
    fgets(fid);
    C{n} = textscan(fid,'%s %s %f %f %f %*s %f %f %f');
    fclose(fid);
%     
    if useInterpolated
        C{n} = cell2mat(C{n}(3:end));
    else
        idxknown = strncmp('known_vector',C{n}{1},12);
        C{n} = cell2mat(C{n}(3:end));
        C{n} = C{n}(idxknown,:);
    end
    
    if ~isempty(SLICE)
        C{n} = C{n}(C{n}(:,3)==SLICE-1, :);
    end
   
end

