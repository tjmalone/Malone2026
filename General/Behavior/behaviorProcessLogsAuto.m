%% behaviorProcessLogsAuto
% Automaticaly process behavioral logs for all mice and generate a
% powerpoints of the behavioral figures.

clear; close all; clc

p1 = pwd;

% set split parameters
folderPtn = '*24*';
tokenLen = 2;
folderTypes = {'5A','5B','6A'};

% set behavior parameters
startNew = 1;
nRuns = 0;
lickType = 'all';


%% Split virmen logs into subfolders

% perform split
splitLogDirs(folderPtn,tokenLen,folderTypes)


%% Run behavior analysis for all mice

% identify mouse folders
folders = dir(folderPtn);

for ii = 1:length(folders)
    cd([folders(ii).folder '\' folders(ii).name '\VirmenLogs'])
    behaviorProcessLogs(startNew,nRuns,lickType)

    close
end

cd(p1)


%% Generate powerpoint

import mlreportgen.ppt.*

% initialize powerpoint
ppt = Presentation(['dailyBehavior_' lickType 'Lick.pptx']);

% make title slide
titleSlide = add(ppt,'Title Slide');
replace(titleSlide,'Title','Mouse Behavior Log');
replace(titleSlide,'Subtitle',...
    string(datetime('today','Format','MMMM d, yyyy')));
replace(titleSlide,'Footer',' ')

% add behavior slides
for ii = 1:length(folders)
    % initialize picture slide
    pictureSlide = add(ppt,'Blank');

    % add mouse ID as title
    tb = TextBox();
    tb.FontSize = '36pt';
    tb.X = '450px';
    tb.Y = '50px';
    tb.Width = '800px';
    add(tb,folders(ii).name);
    add(pictureSlide,tb);

    % add figure
    plot1 = Picture([folders(ii).folder '\' folders(ii).name...
        '\VirmenLogs\training.tif']);
    plot1.Width = '800px';
    plot1.Height = '550px';
    plot1.X = '250px';
    plot1.Y = '150px';
    add(pictureSlide,plot1);

    % remove footer
    replace(pictureSlide,'Footer',' ')
end

cd(p1)

close(ppt);
rptview(ppt);
