function r = figure_mass_conservation

t_start=0;
t_end=1;
time_interval = 0.05;  
N = (t_end - t_start)/time_interval + 1;
A = ones(N,1);

figure;
time = 0:time_interval:t_end;
M = 5.637030e-02 * A;
plot(time', M, 'r');
xlabel('time');
ylabel('The total mass of \rho^n_h');
title('Mass conservation');
savefig('Mass_conservation_32.fig')
close;    