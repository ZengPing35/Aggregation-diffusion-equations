function [rho, free_energy]=nonlocal_nonlinear_aggregation_diffusion_solver_2D(P_partition, T_partition, t_start, t_end, tau, num_control_volume, measure_control_volume, max_iteration_step, tolerence, time_interval)
%% PDE
% \rho_t = \nabla \cdot \Big[\rho \nabla \big(H'(\rho)+V+W*\rho \big) \Big]  in $\Omega \times [0, T]$, \\
% \rho \nabla \big(H'(\rho)+V+W*\rho \big) \cdot \bm{n} = 0 & on $\Gamma \times [0, T]$, \\
% \rho(\cdot,0)=\rho_0 & in $\Omega$.

% We apply the proposed scheme to solve PDE problem in an L-shaped area with 
% %
% \[
% H(\rho) = \frac{1}{m - 1} \rho^m \quad(m = 6), \qquad V = e^{ 0.6 ( (x - 1)^2 + (y - 1)^2) }, \qquad W = e^{  0.6( (x - 1)^2 + (y - 1)^2) }. 
% \]
% %
% The initial data is $\rho_0 = e^{-30 (x^2 + y^2) } + 0.01 $.

%% annotation
%The problem domain is L-shape.

%% basis_type_rho==P1
matrix_size=[num_control_volume num_control_volume];
vector_size=num_control_volume;
P_basis_rho=P_partition;
T_basis_rho=T_partition;

%% Initialize the iteration in time t = tau.
rho_old_time=get_initial_vector('function_rho_initial', P_partition);

%% Free energy at T = 0
FE = function_free_energy(P_partition, num_control_volume, rho_old_time, measure_control_volume);

%% Total mass at T = 0
rho_T0_total_mass=rho_total_mass(rho_old_time, measure_control_volume, num_control_volume);
fprintf('******The total mass at T=0: %d*****\n',rho_T0_total_mass);

%% Figures (T = tau)
% patch
ii = 0;
figure;
X_rho = P_basis_rho(1,:)';
Y_rho = P_basis_rho(2,:)';
patch('Faces',T_basis_rho','Vertices',[X_rho,Y_rho],'facevertexCdata',rho_old_time,'edgecolor','none','facecolor','interp');
title('\rho_h');
colorbar;
colormap(jet);
picturename = strcat('rho_h_patch_',num2str(ii),'.fig');
saveas(gca,picturename,'fig');
close;

% fill3
[m,n] = size(T_basis_rho);
XX = zeros(m,n);
YY = zeros(m,n);
ZZ = zeros(m,n);
for i = 1 : n
    XX(:, i) = P_basis_rho(1, T_basis_rho(:, i));
    YY(:, i) = P_basis_rho(2, T_basis_rho(:, i));
    ZZ(:, i) = rho_old_time(T_basis_rho(:, i));
end
figure;
p = patch;
p.FaceColor = 'interp';
C = ZZ;
fill3(XX, YY, ZZ, C);
title('\rho_h');
colorbar;
colormap(jet);
picturename = strcat('rho_h_fill3_',num2str(ii),'.fig');
saveas(gca,picturename,'fig');
close;

%% Iteration in time.
ii = 1;
nn = 2;
num_time = ( (t_end - t_start) / time_interval ) + 1;
free_energy = zeros(num_time, 1);
free_energy(1) = FE;        %%初始时刻的自由能
N=(t_end-t_start)/tau;
for n=0:N-1           
    current_time=t_start+tau*(n+1);
%     fprintf('*********************Time%d*********************\n',current_time);
    
    %% Assemble the load vector of rho.
    B=assemble_vector_2D_time(rho_old_time, vector_size, num_control_volume, measure_control_volume, tau);
    
    %% Iteration method
     %%Initialize the iteration
    rho_old_n=rho_old_time;    
    for l=1:max_iteration_step
%         fprintf('*********************迭代%d次*********************\n',l);
        
        %% Assemble the matrix A.   
        A=assemble_matrix_nonlinear_2D(P_partition, T_partition, tau, matrix_size, num_control_volume, measure_control_volume, rho_old_n, rho_old_time);
        rho_old=A\B;        
        err_iteration = norm(rho_old_n-rho_old);
%         fprintf('误差%d\n',err_iteration);
        if (err_iteration)<tolerence    
%             fprintf('********迭代%d收敛*******\n',l);
            break;
        end
        if l == max_iteration_step
            fprintf('Not convergence!!!!!!!!!\n');
        end
        rho_old_n=rho_old;
    end
    rho_old_time=rho_old;  
    
    if rem(current_time,time_interval)==0
        %% Total mass at current_time
        rho_T_total_mass=rho_total_mass(rho_old_time, measure_control_volume, num_control_volume);
        fprintf('******The total mass at T=%d: %d*****\n',current_time,rho_T_total_mass);
        
        %% Free energy
        FE = function_free_energy(P_partition, num_control_volume, rho_old_time, measure_control_volume);
        free_energy(nn) = FE;
        nn = nn + 1;
        
        %% save solution
        filename = strcat('rho_old_time_',num2str(ii),'.mat');
        save(filename, 'rho_old_time');
             
        %% Figures----》 numerical solution
        % patch
        figure;
        patch('Faces',T_basis_rho','Vertices',[X_rho,Y_rho],'facevertexCdata',rho_old_time,'edgecolor','none','facecolor','interp');
        title('\rho_h');
        colorbar;
        colormap(jet);
        picturename = strcat('rho_h_patch_',num2str(ii),'.fig');
        saveas(gca,picturename,'fig');
        close;          

        % fill3
        [m,l] = size(T_basis_rho);
        XX = zeros(m,l);
        YY = zeros(m,l);
        ZZ = zeros(m,l);
        for i = 1 : l
            XX(:, i) = P_basis_rho(1, T_basis_rho(:, i));
            YY(:, i) = P_basis_rho(2, T_basis_rho(:, i));
            ZZ(:, i) = rho_old_time(T_basis_rho(:, i));
        end
        figure;
        p = patch;
        p.FaceColor = 'interp';
        C = ZZ;
        fill3(XX, YY, ZZ, C);
        title('\rho_h');
        colorbar;
        colormap(jet);
        picturename = strcat('rho_h_fill3_',num2str(ii),'.fig');
        saveas(gca,picturename,'fig');
        close; 
        %%
        ii=ii+1;
    end
end
rho=rho_old_time;

