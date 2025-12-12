%% statsFull_manual
% Used to generate detailed statistics information for compliance
% with journal manuscript requirements. Can run two-way ANOVA (ordinary or
% repeated measure), Pearson Correlation, and Student's t-test. Functions
% have multiple inputs allowing for flexibility in terms of the input data
% format, paired/unpaired tests, correction for multiple comparisons, and
% repeated measures details. All tests are two-tailed, so individual
% functions must be edited to perform one-tailed analysis. Test name and
% units can be defined.
% 
% Functions:
%   anovaEffectSize - runs ordinary or repeated measure two-way ANOVA. Does
%       not run multiple comparisons. Currently only the group difference
%       is output. Input data is assumed to be data x time and can be two
%       numerical arrays (defaults to repeated measure across time) or two
%       cell arrays (defaults to ordinary). isDep input defines muleiple
%       comparison variables. Requires mes2way.m.
%   
%   corrEffectSize - runs a Pearson correlation between two numerical
%       arrays (if data2 is not empty) or between the mean or a matrix and
%       the indices of the matrix. Input matrix is assumed to be repeat x
%       time. To check that input array is the correct orientation, check n
%       output.
%
%   ttestEffectSize - runs a Student's t-test between two numerical arrays
%       or between the corresponding copumns of two cell arrays. Multiples
%       tests can only be performed with a cell array. All tests will be
%       output in a single row with comma or semicolon delimiters between
%       tests. Can be set to paired/unpaired, corrected/uncorrected for
%       multiple comparisons using the Bonferroni-Holm method, set to only
%       show the results where p is significant, or to perform a ranksum
%       test. A paired test will cause an error if the number of elements
%       is not equal.The nonparametric ranksum test has fewer allowable
%       inputs and provides fewer outputs. Multiple correction requires
%       bonferroni_holm.m.
%
%   ttestAllPairsEffectSize - runs Student's t-tests for pairwise
%       combinations of columns within a single numerical array or cell
%       array. Can be set to paired/unpaired and corrected/uncorrected for
%       multiple comparisons using the Bonferroni-Holm method. A paired
%       test will cause an error if the numebr of elements is not equal.
%       Each test will be output to a different row. The set of pairwise
%       comparisons to perform can be specified (default is all
%       comparisons).Multiple correction requires bonferroni_holm.m.
%
% How to run:
%   1. Run first section to set data1 and data2 as empty numerical or cell
%       arrays.
%   2. Paste data into data1 and data2 as relevant. NaN values and blank
%       values (in cell arrays) are generally ignored or removed in tests.
%   3. In the desired test section, set the exact test name, units for n
%       values, and any relevant input parameters.
%   4. Run the desired test section. Statistics are output in outStats
%
%%%%
% NOTE: This code is designed for flexibility to accommodate a variety of
%   test designs. As such, it is important to check that the data format
%   and all input paramters are set as necessary for the desired test. If
%   possible, it is best to check the output p value against a known
%   correct test for validation.
%%%%
%

clear; clc

% data1 = [];
% data2 = [];

data1 = {};
data2 = {};


%% Perform ANOVA

testName = 'two-way repeated measures ANOVA';

% nUnits = 'cells';
nUnits = 'cells from 11,4 mice';
% nUnits = 'FOV from 11,4 mice';
% nUnits = 'mice';
% nUnits = 'cells from 4,2 mice';
% nUnits = 'sessions from 5 mice';
% nUnits = 'FOV from 4 mice';


outStats = anovaEffectSize(data1,data2,testName,nUnits);
% outStats = anovaEffectSize(data1,data2,testName,nUnits,[1 1]);



%% Perform Pearson Correlation

testName = 'two-tailed Pearson linear correlation';

% nUnits = 'time points';
% nUnits = 'time points from 15 mice';
% nUnits = 'time points from 11 mice';
% nUnits = 'time points from 660 cells from 11 mice';
% nUnits = 'time points from 366 cells from 11 mice';
% nUnits = 'time points from 35 cells from 11 mice';
% nUnits = 'time points from 21 FOV from 11 mice';
% nUnits = 'time points from 8 FOV from 4 mice';

% nUnits = 'time points from 4 mice';
% nUnits = 'time points from 23 cells from 4 mice';
% nUnits = 'time points from 527 cells from 4 mice';
% nUnits = 'time points from 284 cells from 4 mice';
% nUnits = 'time points from 23 cells from 4 mice';
% nUnits = 'time points from 9 FOV from 4 mice';

% nUnits = 'sessions from 15 mice';

% nUnits = 'mice';
nUnits = 'cell pairs';
% nUnits = 'sessions from 6 mice';


outStats = corrEffectSize(data1,data2,testName,nUnits);
% outStats = corrEffectSize(data1',data2',testName,nUnits);


%% Perform Student's t-test

% testName = 'two-tailed unpaired Students t-test';
testName = 'two-tailed unpaired Students t-test, no correction for multiple comparisons';
% testName = 'two-tailed paired Students t-test';
% testName = 'two-tailed paired Students t-test, no correction for multiple comparisons';
% testName = 'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
% testName = 'two-tailed paired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
% testName = 'two-sided Mann-Whitney U test, no correction for multiple comparisons';
% testName = 'one-tailed unpaired Students t-test, no correction for multiple comparisons';


% nUnits = 'cells';
% nUnits = 'cells from 11 mice';
% nUnits = 'cells from 4 mice';
% nUnits = 'cells from 11,4 mice';
% nUnits = 'FOV from 11,4 mice';
% nUnits = 'FOV from 4 mice';
% nUnits = 'mice';
nUnits = 'bins';
% nUnits = 'trials';
% nUnits = 'cell pairs from 11 mice';
% nUnits = 'cell pairs from 4 mice';
% nUnits = 'cells from 6,4 mice';
% nUnits = 'slices from 16,10 mice';
% nUnits = 'sessions from 8,5 mice';


% pair = 1;
pair = 0;

MC = 0;

limitP = 0;
% limitP = 1;

outStats = ttestEffectSize(data1,data2,testName,nUnits,pair,MC,limitP);
% outStats = ttestEffectSize(data1',data2',testName,nUnits,pair,MC,limitP);
% outStats = ttestEffectSize(data1,data2,testName,nUnits,pair,MC,limitP,'ranksum');


%% Perform pairwise Student's t-test

% testName = 'two-tailed unpaired Students t-test';
testName = 'two-tailed unpaired Students t-test, no correction for multiple comparisons';
% testName = 'two-tailed paired Students t-test';
% testName = 'two-tailed paired Students t-test, no correction for multiple comparisons';
% testName = 'two-tailed unpaired Students t-test, Bonferroni-Holm correction for mutiple comparsions';
% testName = 'two-tailed paired Students t-test, Bonferroni-Holm correction for mutiple comparsions';

% nUnits = 'cells';
% nUnits = 'cells from 11 mice';
% nUnits = 'cells from 4 mice';
% nUnits = 'mice';
% nUnits = 'sessions from 8 mice';
% nUnits = 'sessions from 5 mice';
% nUnits = 'FOV from 11,4 mice';
% nUnits = 'FOV from 15 mice';
nUnits = 'bins from 10 sessions';


% pair = 1;
pair = 0;

% MC = 1;
MC = 0;

testIdxs = [1 2;1 3;2 3];
% testIdxs = [1 4;2 5;3 6;2 3;5 6];
% testIdxs = [1 2;3 4;5 6;7 8];
% testIdxs = [];

outStats = ttestAllPairsEffectSize(data1,testName,nUnits,pair,MC,testIdxs);

