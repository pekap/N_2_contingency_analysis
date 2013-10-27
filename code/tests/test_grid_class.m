% test Grid_class.m work 
% you can load any MATPOWER case.
%     for example  loadcase('case2737sop');

clear all;
runcase = loadcase('case300'); % load MATPOWER unreduced Polish case
grid = Grid_class(runcase,0,'case300');

[number_of_0_violations, margin_0_absolute, margin_0_relative,top_0] = grid.N_0_analysis();
grid.N_1_analysis(); % runs N-1 analysis
grid.N_2_analysis('fast'); % runs fast N-2 contingency analysis
grid.N_2_analysis('bruteforce'); % runs bruteforce N-2 contingency analysis

