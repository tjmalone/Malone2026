function imageLapIdx = findImageLaps(folder,svON)
%%

p = pwd;
if nargin>0 && ~isempty(folder)
    cd(folder)
end


%%

load('md.mat','md')

y = md(:,3)*100;
lap = md(:,4)*100;

% find y start and stop
startY = round(y(1));

% find lap start and stop
startLap = round(min(lap))+1;
stopLap = round(max(lap));

load('RunByRun_sig\dfofM_sig.mat','dfofM_sig')
nLaps = size(dfofM_sig{1},1);
if abs(stopLap+1-startLap+1-nLaps)>1 && nLaps~=20
    disp(['Alternate: ' pwd])
    
    load('abf.mat','abf')
    
    y = md(abf.imageIndex,3)*100;
    lap = md(abf.imageIndex,4)*100;
    
    % automatically remove first lap (may be changed) (obsolete)
    % startY = round(y(1));
    startY = 1;
    
    % find lap start and stop
    startLap = round(min(lap))+1;
    stopLap = round(max(lap));
end

lapRuns = startLap:stopLap;
lapIdx = 1:length(lapRuns);

% always remove first lap (may be changed)
% if startY~=0 || startLap==1
    lapRuns = lapRuns(2:end);
    lapIdx = lapIdx(2:end);
% end


% set imageLaps
imageLapIdx = struct();
imageLapIdx.lapIdx = lapIdx;
imageLapIdx.lapRun = lapRuns;


%%

if nargin>1 && svON==1
    save('imageLapIdx.mat','imageLapIdx')
end

cd(p)

end