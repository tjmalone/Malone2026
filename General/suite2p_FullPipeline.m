%% Full Suite2p Pipeline
%
%
%% Motion Correction:
% AnalysisCode\Pipelines\suite2p
%
%   runPipeAuto_dual - Automatically runs motion correction pipeline for
%       all FOVs (locs) for a given mouse.
%   runPipeline_dual - Runs motion correction for a given FOV. Called by
%       runPipeAuto_dual.
%   indivTiffsToStack_dual - Support function that converts individual
%       tifs into stacks.
%   correctMCArea - Support function that corrects all motion corrected
%       files to the smallest common area.
%   correctAbf - calculates additional behavioral paramters from voltage
%       recording if they were skipped initially.
%
%
%% ROI Segmentation/Activity Extraction:
% AnalysisCode\Pipelines\suite2p
%
%   suite2p - performs ROI segmentation and performs initial activity
%       extraction (F). Run in python via Anaconda Prompt. Installation
%       instructions are found in:
%       Z:\commonCodes\Segmentation\Suite2p\updatedInstallInstructions_suite2p
%   suite2pAuto_PostProcess - Automatically runs suite2p post-processing
%       for all suite2p folders at the desired level. This includes
%       determining final ROIs and caluclating of dfof.
%   suite2p_PostProcess - Runs suite2p post-processing for a given suite2p
%       folder. Called by suite2pAuto_PostProcess.
%   dfofCheck - Allows user to manually check ROIs and spatial dfof to
%       identify any irregularities.
%
%
%% Data Backup:
% AnalysisCode\Copying
%
%   robocopy - Most efficient method for copying raw TSeieres data, but can
%       be adapted to copy various file types. Runs from the command
%       terminal, but could potentially be adapted to run in Matlab.
%   copyAllBackupFiles - Automatically copies motion correction and other
%       files to desired drives. Multiple correction, behavior, and run
%       files are copied to a external drive by default. Suite2p folder is
%       copied to the server by default.
%   deleteFiles - Automatically deletes specific files or folders based on
%       inputs. Must be used very carefully to avoid accidental deletions.
%
%
%% Post-Analysis:
% AnalysisCode\PostAnalysis
%
%   copyCues - Copies cue templates to suite2p folders based on use inputs.
%   postAnalysisAuto - Automatically performs post-analysis for all suite2p
%       folders at the desired level. See below for details.
%   postAnalysis_TM - Performs post-analysis for a given suite2p folder.
%       Called by postAnalysisAuto.
%   dfofMCorrelationFull - Calculates spatial dfof, including interpolated
%       and averagesmooth, and calculates run-by-run correlations.
%   pValueRunByRun_dfof_oldSig_noRBR_par - Identifies grid cells using dfof
%   pValueRunByRun_sig_noRBR_par - Identifies grid cells using dfof_sig
%   spatialInfo_dual - Calculates spatial information score
%   speedScoreCalculation - Identifies speed cells.
%   cueCellsAnalysis_dfof_DN_TM - Calculates cue score using dfof.
%   cueCellsAnalysis_sig_DN_TM - Calculates cue score using dfof_sig.
%   postAnalysisCheck - Checks whether all post-analysis ran correctly for
%       all suite2p folders.
%
% Note: alternative versions of some codes allowing for nans in
%   F/dfof/dfof_sig can be found in:
%   AnalysisCode\PostAnalysis\postAnalysis_fix
%
%
%% Behavior Preparation:
% AD_project\Behavior\code and AnalysisCode\Behavior
%
%   trainingLaps - Generates a time course of training data for a mouse
%       over multiple days. Combines sessions of the same type from the
%       same day. After verification, user must copy output to combined
%       excel sheet and then to a matlab structure to generate combined
%       experience data.
%   behaviorProcessLogs - Generates a time course of behavior data for a
%       mouse over multiple days. Is used to determine optimal Day 0 and
%       pre-imaging training days.
% 	defineFiles_behavior - Allows user to define the various sub-categories
% 	    of virmen sessions for behavior analysis. Defines files in 3 stages:
% 	    direct file paths, environment or behavior category by file, robust
% 	    set of indices for various categories.
%
%
%% Alignment:
% AnalysisCode\Alignment and AD_Project\imagingData\code
%
%   defineFilesBasic - Allows user to set genotype and imaging days to be
%       used in alignment. Only one alignment category
%   defineFilesSubtype - Allows user to set genotype and imaging days to be
%       used in alignment. Allows definition of various sub-groups
%   pipeline_align - Contains general outline of all alignment steps. Can be
%       used to run final alignment steps (after manual alignment).
%   alignAutoSuite2p - Performs automated alignment for all FOV. Is set up
%       for the multi-type output of defineFiles, but can be simplified.
%   alignmentSuite2p - Performs automated alignment for a paricular FOV.
%   commonCellIdentification - Performs automated alignment for a
%       particular FOV pair.
%   commonCellIdentificationManual - Performs alignment for a particular
%       FOV pair based on the results of manual alignment.
%   commonCellIdentificationNR - Performs automated alignment for a
%       particular FOV pair using non-rigid alignment.
%   alignCheck - Allows user to define good or bad alignments to determine
%       which pairs of FOV must be manually aligned.
%   runManAlignment - Allows the user to perform manual alignment for a
%       particular FOV pair. Each section should be run separately.
%       Detailed use instructions can be found in
%       "AnalysisCode\Alignment\Manual Alignment Instructions.docx".
%       Requires FOV_manual_alignment_new and alignmentSingle_suite2p.
%   FOV_manual_alignment_new - Prepares manual alignment information for
%       use in the automated alignment code.
%   alignmentSingle_suite2p - Performs alignment for a single FOV pair.
%       Runs using manual alignment (if manual alignment output exists) or
%       uses non-rigid alignment.
%   alignmentPost - Generates a results file combining the alignments for a
%       defined set of alignment FOV. All pairwise alignments must be
%       complete.
%   uniqueAligns - Identifies potential alignments of the desired FOV
%       coverage from an alignment results file.
%   peakAligns - Automatically determines optimal alignments based on a
%       confidence score. Should not be run with alignmentPostCheck.
%   alignmentPostCheck - Allows the user to manually define good and bad
%       alignments. Check must be partly performed in Excel and copied back
%       to Matlab. Should not be run with peakAligns.
%   saveAligns - Removes rejected alignments from input and saves final
%       alignments.
%   genRecallAligns - Generate the alignment results and final alignments
%       for the FE recall.
%
%   Use specific codes:
%   align11Auto/align12Auto: Automates alignment relative to a specific
%       session in order to add extra alignmets to a previously completed
%       alignment analysis
%   alignmentByIdx: Performs alignment for all sessions relative to a
%       specific reference session
%
%
%% Behavior/Activity Preprocessing:
% AnalysisCode\Copying, AnalysisCode\PostAnalysis,
% AnalysisCode\PostAnalysis\Cues, and AD_Project\imagingData\code
%
%   copyM - Copy m.mat into the desired location outside TSeries folder.
%       For use after backup and deletion.
%   moveM/moveMD - Move m.mat or md.mat files into the desired locations
%       (outside TSeries folder or within suite2p folder, respectively).
%       Best run after backup.
%   findImageLapsAuto - Automatically syncronizes imaging lap number with
%       behavior lap number for all sessions from all mice. Calls
%       findImageLaps.
%   findImageLaps - Syncronizes imaging lap number with behavior lap
%       number for a particular imaging session.
%   dfofMCorrelationAll - Calculate dfof and spatial RBR consistency for
%       all sessions from all mice. Reruns 5cm calculation, ensuring
%       consistent code versions and runs at 2.5cm for drift analysis.
%       Calls dfofMCorrelationFull (see above), dfofMCorrelation_clean, and
%       dfofM_RBRbyLoc.
%   dfofM_RBRbyLoc - Calculates sptail RBR consistency for a given session.
%   calculateRBRbyLoc - Alternative code to run dfofM_RBRbyLoc, but without
%       recalculated dfofM.
%   defineFoldersProccessed - Defines suite2p folders based on alignment
%       folders for learning days. Also saves common cell alignment, and
%       defines which session belongs to which group (WT vs. AD) or sex
%       (female vs. male) for imaging analyses.
%   appendTrueDays - Adds a variable define the true day corresponding to
%       each imaging session to the learning folder file.
%       folders for learning days.
%   defineFoldersProccessedRecall - Defines suite2p folders based on
%       alignment folders for recall days. Also saves common cell
%       alignment, and defines which session belongs to which group (WT vs.
%       AD) or sex (female vs. male) for recall imaging analyses.
%	defineMice - Saves indices defining which mice belong to which group (WT
%	    vs. AD) or sex (female vs. male) for behavior analyses.
%   appendTrueDays - Defines "true" match between imaging and behavior
%       sessions, as some mice are missing imaging on some days.
%   uniThreshTM - Defines combined cue score threshold for learning days.
%   uniThreshTMRecall - Defines combined cue score threshold for FE recall
%       days.
%   uniThreshSpeed_TM - Defines combined speed score threshold for learning
%       days.
%   identifySuccessVsFailRuns - Identifies success and fail runs for all
%       sessions for all mice. Calls getLapIdx.
%   getLapIdx - Identifies success and fail runs for a particular session
%   processSuccessVsFailRuns - Performs field identification for success
%       and fail runs and conditional success and fail runs (based on the
%       previous run behavior) for all sessions for all mice. Calls
%       pValueRunByRun_sig_subset.
%   processSuccessN - Performs field identification for the first N success 
%       runs in a novel environment across an arbitrary number of sessions
%       for all mice. Calls pValueRunByRun_sig_subset.
%   processSuccessVsFailJoint - Performs field identification for success
%       and fail runs appended to the next lap for all sessions for all
%       mice. Calls pValueRunByRun_sig_subset.
%   pValueRunByRun_sig_subset - Performs field identification for a
%       particular imaging session using a subset of runs.
%   findGlobalCellTypes - Finds cue, grid, and speed cells among a set of
%       common cells that meet the cell type requirements in a sufficient
%       number of sessions. A copy exists in both the AD and general
%       analysis folders.
%   findGlobalCellTypesRecall - Finds cue, grid, and speed cells among a
%       set of common cells that meet the cell type requirements in
%       reference and recall sessions.
%   findLocalCellTypes - Finds cue, grid, and speed cells among all cells
%       that meet the cell type requirements in individual sessions.
%   RBRFieldIdent_auto - Identifies fields on a run-by-run basis for a
%       specified set of sessions. Is located in and calls several
%       functions within D:\AnalysisCode\Drift
%
% Note: alternative versions of some codes allowing for nans in
%   F/dfof/dfof_sig can be found in:
%   AnalysisCode\PostAnalysis\postAnalysis_fix
%
%

