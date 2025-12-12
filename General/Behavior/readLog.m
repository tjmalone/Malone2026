function logData = readLog(logPath,nColumns)

% open file
fid=fopen(logPath);
% read data
A=fread(fid,inf,'double');
% close file
fclose(fid);


% find number of columns if not specified
if nargin==1
    % assume first item is time stamp, find next number within in 5 seconds
    nColumns = find((A(2:end)>A(1)) & (A(2:end)<(A(1)+5/24/60/60)),1,'first');
end

% reshape data
dataRem = rem(numel(A),nColumns);
if dataRem==0
    logData = reshape(A(1:end),nColumns,[])';
else
    disp('Remainder error')
    logData = reshape(A(1:end-dataRem),nColumns,[])';
end