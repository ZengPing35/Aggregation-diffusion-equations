function b = assemble_vector_2D_time(rho_old, current_time, P_partition, rho_vector_size, rho_num_control_volume, measure_control_volume, tau)
%P_partition: store the coordinates of all the grid points of the partition,not FE.
%T_partition: store the global indices of the grid points of every element for the partition,not FE.
%vertices: the coordinates of all vertices of a triangular element.

%% Initialization
b=zeros(rho_vector_size,1);

%%
for i = 1:rho_num_control_volume    
    vertices_A =  P_partition(:,i); %%当前i控制体的节点坐标
    f=function_f_rho(vertices_A(1), vertices_A(2), current_time);
    b(i) = measure_control_volume(i) / tau * rho_old(i) + f * measure_control_volume(i);
end
