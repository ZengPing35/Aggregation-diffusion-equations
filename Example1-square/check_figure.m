function check_figure
%The problem domain is Circle.
%% Parameters
format short e
basis_type = 1;
basis_index = [1 2 3]; 
Gauss_point_number=9;
max_iteration_step=60;
tolerence=10^(-7);
t_start=0;
t_end=0.2;
time_interval = 0.1;        %%输出图片及计算自由能的时间间隔

%% 保存输出
diary('data.txt');

%% fix tau = 1/80000
tau=1/400;
fprintf('Fix (tau=%d)\n',tau);
        
%% numerical solution (h = 1/4)

% Mesh information for partition
[P_partition,T_partition,Edge] = getMesh_Gmsh('Square_4.msh');
T_partition = T_partition(1:3, :);
num_nodes = size(P_partition,2);            %%三角剖分的节点个数
num_control_volume = num_nodes;             %%控制体的个数

% Exact solution 
rho_exact = get_initial_vector('function_rho_exact', P_partition, t_end);

% Measure of control volume
measure_control_volume = function_measure_control_volume(P_partition, T_partition, num_control_volume);

% numerical solution                                 
[rho, free_energy]=nonlocal_nonlinear_aggregation_diffusion_solver_2D(P_partition, T_partition, t_start, t_end, tau, num_control_volume, measure_control_volume, max_iteration_step, tolerence, time_interval);

% figure of free energy
filename = 'free_energy_4.mat';
save(filename, 'free_energy');

free_energy_x = 0:time_interval:t_end;
plot(free_energy_x', free_energy, 'r');
xlabel('time');
ylabel('free energy');
title('Free energy');
savefig('Free_energy_4.fig')

% L2-norm                 
rho_L2_error_4 = FE_solution_error_triangle_index(rho, 'function_rho_exact', t_end, P_partition, T_partition, basis_type, basis_index, 0, 0, Gauss_point_number)

% H1-norm
rho_H1_error_x = FE_solution_error_triangle_index(rho, 'rho_exact_solution_x_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 1, 0, Gauss_point_number);
rho_H1_error_y = FE_solution_error_triangle_index(rho, 'rho_exact_solution_y_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 0, 1, Gauss_point_number);
rho_H1_error_4=sqrt(rho_H1_error_x^2 + rho_H1_error_y^2)

%% Figures

%% exact solution at T = t_end
% figure;
X_rho = P_partition(1,:)';
Y_rho = P_partition(2,:)';
% patch('Faces',T_partition','Vertices',[X_rho,Y_rho],'facevertexCdata',rho_exact,'edgecolor','none','facecolor','interp');
% title('\rho');
% colorbar;
% colormap(jet);
% savefig('rho_exact_patch_4_t_end.fig');
% close;
% % scatter diagram 散点图
% figure;
% scatter(X_rho,Y_rho,5,rho_exact);                      
% title('\rho');
% colorbar;
% colormap(jet);
% savefig('rho_exact_scatter_4_t_end.fig');
% close;
% pseudo-color image 伪彩色图
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho_exact,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho');
colorbar;
colormap(jet);
savefig('rho_exact_pcolor_4_t_end.fig');
close;
% % contour map 等高线图
% figure;
% contourf(X,Y,Z);      
% title('\rho');
% colorbar;
% colormap(jet);
% savefig('rho_exact_contourf_4_t_end.fig');
% close;
% 3D surface graph 三维曲面图
figure;
surf(X,Y,Z);   
title('\rho_h');
colorbar;
colormap(jet);
savefig('rho_exact_surf_4_t_end.fig');
close;


%% nemerical solution at T = t_end
% figure;
% patch('Faces',T_partition','Vertices',[X_rho,Y_rho],'facevertexCdata',rho,'edgecolor','none','facecolor','interp');
% title('\rho_h');
% colorbar;
% colormap(jet);
% savefig('rho_numerical_patch_4_t_end.fig');
% close;
% % scatter diagram 散点图
% figure;
% scatter(X_rho,Y_rho,5,rho);                      
% title('\rho_h');
% colorbar;
% colormap(jet);
% savefig('rho_numerical_scatter_4_t_end.fig');
% close;
% pseudo-color image 伪彩色图
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho_h');
colorbar;
colormap(jet);
savefig('rho_numerical_pcolor_4_t_end.fig');
close;
% % contour map 等高线图
% figure;
% contourf(X,Y,Z);      
% title('\rho_h');
% colorbar;
% colormap(jet);
% savefig('rho_numerical_contourf_4_t_end.fig');
% close;
% 3D surface graph 三维曲面图
figure;
surf(X,Y,Z);   
title('\rho_h');
colorbar;
colormap(jet);
savefig('rho_numerical_surf_4_t_end.fig');
close;
