%% saveFilteredRoisReelCal
% saves a new zip file with imageJ rois.
% provide the indeces of the rois you want in the new zip file.
% 
% Main use: save new rois with dim and small cells removed


% for this version
    % input rois:
        % reelin_neurodegen_brightCal_v3.zip
        % reelin_neurodegen_brightCal.zip
    % output rois:
        % reelin_neurodegen_brightCal_v4

clear; close all; clc

%%%%% set values based on cutoffs
foldVersion = 'Z:\SNM\labMembers\KC\AD_Project\Histology\reelCalNeurodegen\FinalVersion\remBorderCells';
finalNameEnding = '_neurodegen_brightCal_v4';
%%%%%

types = {'reelin','calbindin'};

% load indeces
% cd([foldVersion '/' foldVersName])
cd(foldVersion)
load('removeBorderCells.mat')

fNum = length(folds);
% plotTestIndex = [2 7 8 9]; % mice to plot
plotTestIndex = [1:fNum]; % mice to plot

% folder where code is saved
baseFolder = 'D:\MATLAB\Code\AnalysisCode\Analysis\HistologyCode\ReelCalNeurodegen\withRemoveBorderCells';

% copy self to current directory
copyfile([baseFolder '\saveFilteredRoisReelCal.m'],'saveFilteredRoisReelCal.m');

% temp fold location
foldTemp = pwd;

%%
 for cc = 1:length(types) % number of cell types
    zipName = [types{cc} finalNameEnding '.zip'];
    
    for f = 1:fNum
    % for f = 2

        if ismember(f,plotTestIndex)
            % Define paths
            inputZip = [folds{f} '\' cells{f,cc}];
            outputZip = [folds{f} '\' zipName];
            
            % Define the indices of the ROIs you want to keep (1-based index)
            selectedIndices = cellsInsideBorder{f,cc};  % Modify as needed
            
            % - get zp file contents
           cstrFilenames_short = listzipcontents_rois(inputZip);
           
           % - Unzip the file into a temporary directory
           unzip(inputZip, 'temp_rois');
           
           for (nFileIndex = 1:length(cstrFilenames_short))
              cstrFilenames{nFileIndex,1} = ['temp_rois/' char(cstrFilenames_short(nFileIndex, 1))];
              newName{nFileIndex,1} = ['temp_rois/' sprintf('%05d', nFileIndex) '.roi'];
              movefile(cstrFilenames{nFileIndex,1}, newName{nFileIndex,1})
           end
            
            % List all extracted ROI files
            roiFiles = dir(fullfile('temp_rois', '*.roi'));
            roiNames = {roiFiles.name};  % Store file names
            
            % Select the desired ROIs
            selectedRois = roiNames(selectedIndices);
            
            % Create a new zip file with only the selected ROIs
            % zip(outputZip, fullfile('temp_rois', selectedRois));
            cd('temp_rois')
            zip(outputZip, selectedRois);
            cd(foldTemp)
            
            % Clean up (optional)
            rmdir('temp_rois', 's');  % Delete the temporary folder and its contents
            
            disp(['Saved selected ROIs to: ', outputZip]);
        end
    end
 end




 function [filelist] = listzipcontents_rois(zipFilename)
      
      % listzipcontents_rois - FUNCTION Read the file names in a zip file
      %
      % Usage: [filelist] = listzipcontents_rois(zipFilename)
      
      % - Import java libraries
      import java.util.zip.*;
      import java.io.*;
      
      % - Read file list via JAVA object
      filelist={};
      in = ZipInputStream(FileInputStream(zipFilename));
      entry = in.getNextEntry();
      
      % - Filter ROI files
      while (entry~=0)
         name = entry.getName;
         if (name.endsWith('.roi'))
            filelist = cat(1,filelist,char(name));
         end;
         entry = in.getNextEntry();
      end;
      
      % - Close zip file
      in.close();
   end