function [rho, free_energy]=nonlocal_nonlinear_aggregation_diffusion_solver_2D(P_partition, T_partition, t_start, t_end, tau, num_control_volume, measure_control_volume, max_iteration_step, tolerence, time_interval)
%% PDE
% \rho_t = \nabla \cdot \Big[\rho \nabla \big(H'(\rho)+V+W*\rho \big) \Big]  in $\Omega \times [0, T]$, \\
% \rho \nabla \big(H'(\rho)+V+W*\rho \big) \cdot \bm{n} = 0 & on $\Gamma \times [0, T]$, \\
% \rho(\cdot,0)=\rho_0 & in $\Omega$.

%% annotation
%The problem domain is Circle.

%% basis_type_rho==P1
matrix_size=[num_control_volume num_control_volume];
vector_size=num_control_volume;
P_basis_rho=P_partition;
T_basis_rho=T_partition;

%% Initialize the iteration in time t = tau.
rho_old_time=get_initial_vector('function_rho_exact', P_partition, 0);

%% Free energy at T = 0
FE = function_free_energy(P_partition, num_control_volume, rho_old_time, measure_control_volume);
fprintf('******当前时刻T=0的自由能%d*****\n', FE);

%% Total mass at T = 0
% rho_T0_total_mass=rho_total_mass(rho_old_time, measure_control_volume, num_control_volume);
% fprintf('******初始时刻n=0 n的总质量%d*****\n',rho_T0_total_mass);

%% Figures (T = tau)
% ii = 0;
% figure;
% X_rho = P_basis_rho(1,:)';
% Y_rho = P_basis_rho(2,:)';
% patch('Faces',T_basis_rho','Vertices',[X_rho,Y_rho],'facevertexCdata',rho_old_time,'edgecolor','none','facecolor','interp');
% title('\rho_h');
% colorbar;
% colormap(jet);
% picturename = strcat('rho_h_',num2str(ii),'.fig');
% saveas(gca,picturename,'fig');
% close;

%% Iteration in time.
ii = 1;
nn = 2;
num_time = int16(( (t_end - t_start) / time_interval ) + 1);
free_energy = zeros(num_time, 1);
free_energy(1) = FE;        %%初始时刻的自由能
N=(t_end-t_start)/tau;
for n=0:N-1           
    current_time=t_start+tau*(n+1);
%     fprintf('*********************时间%d*********************\n',current_time);
    
    %% Assemble the load vector of rho. 
    B=assemble_vector_2D_time(rho_old_time, current_time, P_partition, vector_size, num_control_volume, measure_control_volume, tau);
    
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
%         rho_T_total_mass=rho_total_mass(rho_old, measure_control_volume, num_control_volume);
%         fprintf('******当前时刻T=%d 第%d次迭代的总质量%d*****\n',current_time, l, rho_T_total_mass);
%         if l == 1
%             rho_T0_total_mass=rho_total_mass(rho_old, measure_control_volume, num_control_volume);
%             fprintf('******当前时刻T=%d 第%d次迭代的总质量%d*****\n',current_time, l, rho_T0_total_mass);
%         end
        if (err_iteration)<tolerence    
%             fprintf('********迭代%d收敛*******\n',l);
%             rho_T_total_mass=rho_total_mass(rho_old, measure_control_volume, num_control_volume);
%             fprintf('******当前时刻T=%d 第%d次迭代的总质量%d*****\n',current_time, l, rho_T_total_mass);
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
%         rho_T_total_mass=rho_total_mass(rho_old_time, measure_control_volume, num_control_volume);
%         fprintf('******当前时刻T=%d n的总质量%d*****\n',current_time,rho_T_total_mass);
        
        %% Free energy
        FE = function_free_energy(P_partition, num_control_volume, rho_old_time, measure_control_volume);
        free_energy(nn) = FE;
        fprintf('******当前时刻T=%d的自由能%d*****\n',current_time, FE);
        nn = nn + 1;
        
        %% save solution
%         filename = strcat('rho_old_time_',num2str(ii),'.mat');
%         save(filename, 'rho_old_time');
        
        %% Figures 
%         %% exact solution
% %         % patch
%         rho_exact=get_initial_vector('function_rho_exact', P_partition, current_time);
% %         figure;
%         X_rho = P_basis_rho(1,:)';
%         Y_rho = P_basis_rho(2,:)';
% %         patch('Faces',T_basis_rho','Vertices',[X_rho,Y_rho],'facevertexCdata',rho_exact,'edgecolor','none','facecolor','interp');
% %         title('\rho');
% %         colorbar;
% %         colormap(jet);
% %         picturename = strcat('rho_patch_',num2str(ii),'.fig');
% %         saveas(gca,picturename,'fig');
% %         close;
% %         % pseudo-color image 伪彩色图
% %         figure;
%         [XX,YY,ZZ]=griddata(X_rho,Y_rho,rho_exact,linspace(-1,1)',linspace(-1,1),'v4');  %插值
% %         pcolor(XX,YY,ZZ);
% %         shading interp;              
% %         title('\rho');
% %         colorbar;
% %         colormap(jet);
% %         picturename = strcat('rho_pcolor_',num2str(ii),'.fig');
% %         saveas(gca,picturename,'fig');
% %         close;
%        % 3D surface graph 三维曲面图
%         figure;
%         surf(XX,YY,ZZ);   
%         title('\rho');
%         colorbar;
%         colormap(jet);
%         picturename = strcat('rho_surf_',num2str(ii),'.fig');
%         saveas(gca,picturename,'fig');
%         close;
%        
%         %% numerical solution
% %         % patch
% %         figure;
% %         patch('Faces',T_basis_rho','Vertices',[X_rho,Y_rho],'facevertexCdata',rho_old_time,'edgecolor','none','facecolor','interp');
% %         title('\rho_h');
% %         colorbar;
% %         colormap(jet);
% %         picturename = strcat('rho_h_patch_',num2str(ii),'.fig');
% %         saveas(gca,picturename,'fig');
% %         close;    
% %         % pseudo-color image 伪彩色图
% %         figure;
%         [XX,YY,ZZ_h]=griddata(X_rho,Y_rho,rho_old_time,linspace(-1,1)',linspace(-1,1),'v4');  %插值
% %         pcolor(XX,YY,ZZ_h);
% %         shading interp;              
% %         title('\rho_h');
% %         colorbar;
% %         colormap(jet);
% %         picturename = strcat('rho_h_pcolor_',num2str(ii),'.fig');
% %         saveas(gca,picturename,'fig');
% %         close;
%        % 3D surface graph 三维曲面图
%         figure;
%         surf(XX,YY,ZZ_h);   
%         title('\rho_h');
%         colorbar;
%         colormap(jet);
%         picturename = strcat('rho_h_surf_',num2str(ii),'.fig');
%         saveas(gca,picturename,'fig');
%         close;
%         %%
%         ii=ii+1;
    end
end
rho=rho_old_time;

