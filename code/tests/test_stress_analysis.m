% stress analysis for case300

clear all;
runcase = loadcase('case300'); % load MATPOWER unreduced Polish case
grid = Grid_class(runcase,'case300_stress');
runcase = grid.rnc;
maximum_stress = Stress_class.maximum_stress(runcase);
rnc = Stress_class.stress(runcase,maximum_stress);
grid = Grid_class(rnc,'case300_stress');
grid = grid.N_1_analysis();
grid = grid.N_2_analysis('bruteforce');
rnc.branch(:,6) = grid.lim;
%grid = Grid_class(rnc,'case300_stress');

n = 100; % number of measurements
filtered = zeros(n,1);
for i=1:n
   grid1 = grid;
   coef = i*(maximum_stress - 1) / n + 1;
   grid1.f = grid1.f * coef/maximum_stress;
   grid1 = grid1.N_1_analysis();
   grid1 = grid1.N_2_analysis('fast');
   filtered(i) = grid1.filtered_size;
   final(i) = Sz.r(grid1.brute_cont);
   time_fast(i) = grid1.t_fast;
   time_brute(i) = grid1.t_brute;
   %grid = grid.N_2_analysis('bruteforce');
   %time_brute2(i) = grid.t_brute;
end


plot(1:n,filtered)
hold on;
plot(1:n,final)
plot(1:n,time_fast*1000,'r')
plot(1:n,(time_fast-time_brute)*10000,'r')