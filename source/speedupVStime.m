% speedup
clear all; clc; close all;

% #threads
threads = [1 2 3 4 5 6 7 8];

% Experimental data

% poisson; 1 thread; 2 threads; 3 thread; 4 thread; 5 thread; 6 thread; 7, thread, thread 8;
poisson_time = [16.423573 9.232522 6.162515 4.867531 4.037279 3.427179 3.101819 2.758126
                16.832716 9.025833 5.991352 5.032488 4.005653 3.426717 3.091921 2.717910
                16.513080 9.254303 6.031612 5.077150 3.984197 3.426052 3.044845 2.739263];  
poisson_time_mean = mean(poisson_time);

% Calculating speedup, S = Tinit/Tnew
speedup = poisson_time_mean(1)./poisson_time_mean(:);

figure
plot(threads,poisson_time_mean,threads,speedup, 'LineWidth',2)
set(gca,'FontSize',13, 'FontWeight', 'bold');
title(['Time cost and speedup']);
ylabel('Speedup and time cost') % label left y-axis
xlabel('Threads') % label x-axis
legend('Time cost','Speedup')
grid on
