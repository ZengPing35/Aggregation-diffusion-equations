function main
%% The problem domain is Annulus.
% In this example, we explore the progression of density across the annulus. In our numerical test, we take the initial data $\rho_0 = e^{ -30  (x^2 + y^2) } + 0.01$. 
% The computational domain is $\Omega = \big\{ (x, y) | 1 \le x^2 + y^2 \le 4 \big\}$. The mesh size and time step are set as $\tau = \frac{1}{20000}$ and $h = 0.0625$, respectively. We take
% %
% \[
% H(\rho) = \frac{1}{m - 1} \rho^m \quad (m = 2), \qquad V =  e^{x^2 + (y + 1)^2}, \qquad W = e^{x^2 + (y + 1)^2}.
% \]
% %

%% Parameters
format short e
max_iteration_step=60;
tolerence=10^(-6);
t_start=0;
t_end=3;
time_interval = 0.1;        %%���ͼƬ�����������ܵ�ʱ����

%% �������
diary('data.txt');

%% fix tau = 1/20000
tau=1/20000;
fprintf('Fix (tau=%d)\n',tau);
        
%% numerical solution (h = 1/16)

% Mesh information for partition
[P_partition,T_partition,Edge] = getMesh_Gmsh('Annulus_16.msh');
T_partition = T_partition(1:3, :);
num_nodes = size(P_partition,2);            %%�����ʷֵĽڵ����
num_control_volume = num_nodes;             %%������ĸ���

% Measure of control volume
measure_control_volume = function_measure_control_volume(P_partition, T_partition, num_control_volume);

% numerical solution                                 
[rho, free_energy]=nonlocal_nonlinear_aggregation_diffusion_solver_2D(P_partition, T_partition, t_start, t_end, tau, num_control_volume, measure_control_volume, max_iteration_step, tolerence, time_interval);

% figure of free energy
figure;
free_energy_x = 0:time_interval:t_end;
plot(free_energy_x', free_energy, 'r');
xlabel('time');
ylabel('free energy');
title('Free energy');
savefig('Free_energy_16.fig')
close;     

%% Figures
% nemerical solution at T = t_end
figure;
X_rho = P_partition(1,:)';
Y_rho = P_partition(2,:)';
patch('Faces',T_partition','Vertices',[X_rho,Y_rho],'facevertexCdata',rho,'edgecolor','none','facecolor','interp');
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_patch_16_t_end.fig');
close;

% fill3
[m,l] = size(T_partition);
XX = zeros(m,l);
YY = zeros(m,l);
ZZ = zeros(m,l);
for i = 1 : l
    XX(:, i) = P_partition(1, T_partition(:, i));
    YY(:, i) = P_partition(2, T_partition(:, i));
    ZZ(:, i) = rho(T_partition(:, i));
end
figure;
p = patch;
p.FaceColor = 'interp';
C = ZZ;
fill3(XX, YY, ZZ, C);
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_fill3_16_t_end.fig');
close; 