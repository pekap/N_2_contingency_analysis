% test Grid_class.m work 
clear all;

pol_grid = Grid_class('sop_reduced.mat');
if (pol_grid.err==0)
   pol_grid.dcopf();
else
   fprintf('there were errors "%s" before - cant run OPF',pol_grid.error);
end
