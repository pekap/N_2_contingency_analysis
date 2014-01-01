% stress analysis plot results
load('case300_100_simulations');
close all;
figure();
x=(1:n)*maximum_stress/n;
[AX,H1,H2] = plotyy(x,[filtered,final'],x,[time_fast],'semilogy','semilogy');
set(get(AX(1),'Ylabel'),'String','C_2^{filtered} set sizes');
set(get(AX(2),'Ylabel'),'String','Completion times');
set(get(AX(1),'Xlabel'),'String','Load coefficient');
set(AX(1),'ycolor','b');
set(AX(2),'ycolor','r');
set(findall(AX,'type','text'),'fontSize',14,'fontname', 'Roboto');
set(H1(1),'LineStyle','--');
set(H1(2),'LineStyle','-')
set(H1(1),'color','b');
set(H1(2),'color','b');
set(H2,'LineStyle','--')
set(H2,'color','r');
legend('Size of  the set C_2^{filtered}','Size of  the set C_2^{final}','Completion time');


figure();
[AX,H1,H2] = plotyy(x,[filtered'./final],x,[time_fast;time_fast-time_brute],'plot','plot');
set(get(AX(1),'Ylabel'),'String','|C_2^{filtered}|/|C_2^{final}|');
set(get(AX(2),'Ylabel'),'String','Completion time');
set(get(AX(1),'Xlabel'),'String','Load coefficient');
set(AX(1),'ycolor','b');
set(AX(2),'ycolor','r');
set(findall(AX,'type','text'),'fontSize',14);
set(H2(1),'LineStyle','--');
set(H2(2),'LineStyle','-')
set(H2(1),'color','r');
set(H2(2),'color','r');
xlim(AX(1),[6 15])
xlim(AX(2),[6 15])
set(H1,'LineStyle','--')
set(H1,'color','b');
legend('|C_2^{filtered}|/|C_2^{final}|','Completion time of the pruning loop','Completion time of the complete search');


semilogy(1:n,(time_fast-time_brute),'g');
%set(AX,'XTickLabel','1|10|100');
set(findall(AX,'type','text'),'fontSize',12,'fontname', 'Monaco');
% semilogy(1:n,filtered)
% hold on;
% semilogy(1:n,final)
% figure();
% plot(1:n,filtered./final');
% figure();
% semilogy(1:n,time_fast,'r')
% hold on;
% semilogy(1:n,(time_fast-time_brute),'g')