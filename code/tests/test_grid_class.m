% test Grid_class.m work 
% rncreduced
% case2737sop

clear all;
 
pol_reduced = load('sop_reduced_vova'); % load Polish grid reduced by Vova
pol_full = loadcase('case300'); % load MATPOWER unreduced Polish case
% 'case2737sop'
grid = Grid_class(pol_full,0,'case300');
[number_of_0_violations, margin_0_absolute, margin_0_relative,top_0] = grid.N_0_analysis();

%grid.N_1_analysis();
grid.N_2_analysis('fast');
grid.N_2_analysis('bruteforce');
%pol_grid.dcopf();

