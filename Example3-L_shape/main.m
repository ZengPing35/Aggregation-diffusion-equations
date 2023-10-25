function main

%% PDE
% \rho_t = \nabla \cdot \Big[\rho \nabla \big(H'(\rho)+V+W*\rho \big) \Big]  in $\Omega \times [0, T]$, \\
% \rho \nabla \big(H'(\rho)+V+W*\rho \big) \cdot \bm{n} = 0 & on $\Gamma \times [0, T]$, \\
% \rho(\cdot,0)=\rho_0 & in $\Omega$.

%The problem domain is L-shape
% We apply the proposed scheme to solve PDE problem in an L-shaped area with 
% %
% \[
% H(\rho) = \frac{1}{m - 1} \rho^m \quad(m = 6), \qquad V = e^{ 0.6 ( (x - 1)^2 + (y - 1)^2) }, \qquad W = e^{  0.6( (x - 1)^2 + (y - 1)^2) }. 
% \]
% %
% The initial data is $\rho_0 = e^{-30 (x^2 + y^2) } + 0.01 $.

%% Parameters
format short e
max_iteration_step=60;
tolerence=10^(-6);
t_start=0;
t_end=1;
time_interval = 0.05;        %%输出图片及计算自由能的时间间隔

%% 保存输出
diary('data.txt');

%% fix tau = 1/60000
tau=1/60000;
fprintf('Fix (tau=%d)\n',tau);
        
%% numerical solution (h = 1/32)

% Mesh information for partition
[P_partition,T_partition,Edge] = getMesh_Gmsh('L_shape_32.msh');
T_partition = T_partition(1:3, :);
num_nodes = size(P_partition,2);            %%三角剖分的节点个数
num_control_volume = num_nodes;             %%控制体的个数

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
savefig('Free_energy_32.fig')
close;     

%% Figures
% nemerical solution at T = t_end
figure;
X_rho = P_partition(1,:)';
Y_rho = P_partition(2,:)';
patch('Faces',T_partition','Vertices',[X_rho,Y_rho],'facevertexCdata',rho,'edgecolor','none','facecolor','interp');
title('\rho_h');
colorbar;
colormap(jet);
savefig('rho_numerical_patch_32_t_end.fig');
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
colorbar;
colormap(jet);
savefig('rho_numerical_fill3_32_t_end.fig');
close; 