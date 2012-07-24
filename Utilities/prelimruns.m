
clear all
close all
more off
format long g


n_july06=   [1   ,2   ,4   ,6   ,8   ,10 ,12 ,14 ,16 ,18 ];
time_july06=[39.2,22.7,11.6,8.85,6.17,5.6,4.4,4.0,3.5,3.3];
n_nov06=    [1   ,4   ,8   ,12 ,16  ,24  ,32  ,40  , 48 ];
time_nov06= [23.5,7.21,4.35,3.3,2.24,1.81,1.56,1.39,1.31];
time_nov06b=[10.8,3.53,1.84,1.4,1.05,.846,.732,.631,.644];
ideal = 10/10*ones(1,9);


speedup_july06 = time_july06(1)*time_july06.^-1;
speedup_nov06 = time_nov06(1)*time_nov06.^-1;
speedup_nov06b = time_nov06b(1)*time_nov06b.^-1;

figure(1)
plot(n_july06,speedup_july06,'k','linewidth',3)
hold
plot(n_nov06,speedup_nov06,'b','linewidth',3)
plot(n_nov06,speedup_nov06b,'g','linewidth',3)
plot(n_nov06,n_nov06,'r--','linewidth',3)
xlabel('Number of Processors','Fontsize',24)
ylabel('Speedup','Fontsize',24)
set(get(gcf,'CurrentAxes'),'FontSize',24)
legend('(july 2006) 7754 dof 3.0Ghz','(nov 2006) 10161 dof 2.66Ghz','(nov 2006) 10161 dof 2.66Ghz', 'linear','location','northwest')
print -djpeg speedupplotcompare

figure(2)
plot(n_july06,time_july06,'k','linewidth',3)
hold
plot(n_nov06,time_nov06,'b','linewidth',3)
plot(n_nov06,time_nov06b,'g','linewidth',3)
plot(n_nov06,ideal,'r--','linewidth',3)
xlabel('Number of Processors','Fontsize',24)
ylabel('Execution Time (sec)','Fontsize',24)
set(get(gcf,'CurrentAxes'),'FontSize',24)
legend('(july 2006) 7754 dof 3.0Ghz','(nov 2006) 10161 dof 2.66Ghz','(nov 2006) 10161 dof 2.66Ghz','Required','location','northeast')
print -djpeg executiontimecompare

figure(3)
plot(n_nov06,time_nov06b,'k','linewidth',3)
hold
plot(n_nov06,ideal,'r--','linewidth',3)
xlabel('Number of Processors','Fontsize',24)
ylabel('Execution Time (sec)','Fontsize',24)
set(get(gcf,'CurrentAxes'),'FontSize',24)
legend('10161 dof 2.66Ghz','location','northeast')
print -djpeg executiontime
