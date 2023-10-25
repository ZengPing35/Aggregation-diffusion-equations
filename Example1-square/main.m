function main
%% The problem domain is Square.
% We consider the initial value problem with a source term,
% %
% \begin{subnumcases}{\label{eq: ADE_F}}
% \nonumber
% \rho_t - \nabla \cdot \Big[\rho \nabla \big(H'(\rho)+V+W*\rho \big) \Big] = F & in $\Omega \times [0, T]$, \\
% \nonumber
% \rho \nabla \big(H'(\rho)+V+W*\rho \big) \cdot \bm{n} = 0 & on $\partial \Omega \times [0, T]$,\\
% \nonumber
% H(\rho) = \frac{1}{m - 1}\rho^m, \qquad H'(\rho) = \frac{m}{m - 1} \rho^{m - 1}, \\
% \nonumber
% V(x, y) = x^2 y^2 (1-x)^2 (1-y)^2, \\
% \nonumber
% W(x, y) = - 0.001 x^2 y^2,
% \end{subnumcases}
% %
% where we take $m = 3$, $\Omega=(0,1)^2$ and $T=0.2$. The exact solution is given by: 
% %
% \[
% \rho = e^{-t} x^2 y^2 (1-x)^2 (1-y)^2 + 0.1.
% \]
% %
% The source term $F$ is computed by substituting the above exact solution into the equation of $\rho$.

%% Parameters
format short e
basis_type = 1;
basis_index = [1 2 3]; 
Gauss_point_number=9;
max_iteration_step=60;
tolerence=10^(-7);
t_start=0;
t_end=0.2;
time_interval = 0.01;        %%输出图片及计算自由能的时间间隔

%% 保存输出
diary('data_fix_tau.txt')

%% fix tau = 1/10000
tau=1/20000;
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
figure;
free_energy_x = 0:time_interval:t_end;
plot(free_energy_x', free_energy, 'r');
xlabel('time');
ylabel('free energy');
title('Free energy');
savefig('Free_energy_4.fig')
close;

% L2-norm                 
rho_L2_error_FVM_4=FVM_solution_L2error(num_control_volume, measure_control_volume, rho_exact, rho)

% H1-norm
rho_H1_error_x = FE_solution_error_triangle_index(rho, 'rho_exact_solution_x_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 1, 0, Gauss_point_number);
rho_H1_error_y = FE_solution_error_triangle_index(rho, 'rho_exact_solution_y_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 0, 1, Gauss_point_number);
rho_H1_error_4=sqrt(rho_H1_error_x^2 + rho_H1_error_y^2)

%% Figures
%% exact solution at T = t_end
X_rho = P_partition(1,:)';
Y_rho = P_partition(2,:)';
% pseudo-color image 
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho_exact,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_exact_pcolor_4_t_end.fig');
close;
% 3D surface graph 
figure;
surf(X,Y,Z);   
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_exact_surf_4_t_end.fig');
close;

% nemerical solution at T = t_end
% pseudo-color image 
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_pcolor_4_t_end.fig');
close;
% 3D surface graph 
figure;
surf(X,Y,Z);   
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_surf_4_t_end.fig');
close;

%% numerical solution (h = 1/8)

% Mesh information for partition
[P_partition,T_partition,Edge] = getMesh_Gmsh('Square_8.msh');
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
filename = 'free_energy_8.mat';
save(filename, 'free_energy');
figure;
free_energy_x = 0:time_interval:t_end;
plot(free_energy_x', free_energy, 'r');
xlabel('time');
ylabel('free energy');
title('Free energy');
savefig('Free_energy_8.fig')
close;

% L2-norm                 
rho_L2_error_FVM_8=FVM_solution_L2error(num_control_volume, measure_control_volume, rho_exact, rho)
r=(log(rho_L2_error_FVM_4)-log(rho_L2_error_FVM_8))/(log(2));
fprintf('***Convergence order of rho_L2 for h=%f***\n',r)

% H1-norm
rho_H1_error_x = FE_solution_error_triangle_index(rho, 'rho_exact_solution_x_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 1, 0, Gauss_point_number);
rho_H1_error_y = FE_solution_error_triangle_index(rho, 'rho_exact_solution_y_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 0, 1, Gauss_point_number);
rho_H1_error_8=sqrt(rho_H1_error_x^2 + rho_H1_error_y^2)
r=(log(rho_H1_error_4)-log(rho_H1_error_8))/(log(2));
fprintf('***Convergence order of rho_H1 for h=%f***\n',r)

%% Figures
%% exact solution at T = t_end
X_rho = P_partition(1,:)';
Y_rho = P_partition(2,:)';
% pseudo-color image 
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho_exact,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_exact_pcolor_8_t_end.fig');
close;
% 3D surface graph 
figure;
surf(X,Y,Z);   
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_exact_surf_8_t_end.fig');
close;

%% nemerical solution at T = t_end
% pseudo-color image 
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_pcolor_8_t_end.fig');
close;
% 3D surface graph 
figure;
surf(X,Y,Z);   
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_surf_8_t_end.fig');
close;

%% numerical solution (h = 1/16)

% Mesh information for partition
[P_partition,T_partition,Edge] = getMesh_Gmsh('Square_16.msh');
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
filename = 'free_energy_16.mat';
save(filename, 'free_energy');
figure;
free_energy_x = 0:time_interval:t_end;
plot(free_energy_x', free_energy, 'r');
xlabel('time');
ylabel('free energy');
title('Free energy');
savefig('Free_energy_16.fig')
close;

% L2-norm                 
rho_L2_error_FVM_16=FVM_solution_L2error(num_control_volume, measure_control_volume, rho_exact, rho)
r=(log(rho_L2_error_FVM_8)-log(rho_L2_error_FVM_16))/(log(2));
fprintf('***Convergence order of rho_L2 for h=%f***\n',r)

% H1-norm
rho_H1_error_x = FE_solution_error_triangle_index(rho, 'rho_exact_solution_x_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 1, 0, Gauss_point_number);
rho_H1_error_y = FE_solution_error_triangle_index(rho, 'rho_exact_solution_y_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 0, 1, Gauss_point_number);
rho_H1_error_16=sqrt(rho_H1_error_x^2 + rho_H1_error_y^2)
r=(log(rho_H1_error_8)-log(rho_H1_error_16))/(log(2));
fprintf('***Convergence order of rho_H1 for h=%f***\n',r)

%% Figures
%% exact solution at T = t_end
X_rho = P_partition(1,:)';
Y_rho = P_partition(2,:)';
% pseudo-color image 
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho_exact,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_exact_pcolor_16_t_end.fig');
close;
% 3D surface graph 
figure;
surf(X,Y,Z);   
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_exact_surf_16_t_end.fig');
close;

%% nemerical solution at T = t_end
% pseudo-color image 
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_pcolor_16_t_end.fig');
close;
% 3D surface graph 
figure;
surf(X,Y,Z);   
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_surf_16_t_end.fig');
close;


%% numerical solution (h = 1/32)

% Mesh information for partition
[P_partition,T_partition,Edge] = getMesh_Gmsh('Square_32.msh');
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
filename = 'free_energy_32.mat';
save(filename, 'free_energy');
figure;
free_energy_x = 0:time_interval:t_end;
plot(free_energy_x', free_energy, 'r');
xlabel('time');
ylabel('free energy');
title('Free energy');
savefig('Free_energy_32.fig');
close;

% L2-norm                 
rho_L2_error_FVM_32=FVM_solution_L2error(num_control_volume, measure_control_volume, rho_exact, rho)
r=(log(rho_L2_error_FVM_16)-log(rho_L2_error_FVM_32))/(log(2));
fprintf('***Convergence order of rho_L2 for h=%f***\n',r)

% H1-norm
rho_H1_error_x = FE_solution_error_triangle_index(rho, 'rho_exact_solution_x_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 1, 0, Gauss_point_number);
rho_H1_error_y = FE_solution_error_triangle_index(rho, 'rho_exact_solution_y_derivative', t_end, P_partition, T_partition, basis_type, basis_index, 0, 1, Gauss_point_number);
rho_H1_error_32=sqrt(rho_H1_error_x^2 + rho_H1_error_y^2)
r=(log(rho_H1_error_16)-log(rho_H1_error_32))/(log(2));
fprintf('***Convergence order of rho_H1 for h=%f***\n',r)

%% Figures
%% exact solution at T = t_end
X_rho = P_partition(1,:)';
Y_rho = P_partition(2,:)';
% pseudo-color image 
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho_exact,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_exact_pcolor_32_t_end.fig');
close;
% 3D surface graph 
figure;
surf(X,Y,Z);   
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_exact_surf_32_t_end.fig');
close;

%% nemerical solution at T = t_end
% pseudo-color image 
figure;
[X,Y,Z]=griddata(X_rho,Y_rho,rho,linspace(0,1)',linspace(0,1),'v4');  %插值
pcolor(X,Y,Z);
shading interp;              
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_pcolor_32_t_end.fig');
close;
% 3D surface graph 
figure;
surf(X,Y,Z);   
title('\rho_h');
xlabel('x'); ylabel('y');
colorbar;
colormap(jet);
savefig('rho_numerical_surf_32_t_end.fig');
close;

%% 保存输出
diary off;