function r=function_Discrete_convolution_free_energy(P_partition, vertices_i, num_control_volume, measure_control_volume, rho_old_time, i)
% i:第i个控制体标号

r=0;
for n = 1:num_control_volume
    vertices = P_partition(:,n);    
    r = r + (1/2) * measure_control_volume(n) * function_W((vertices_i(1,1) - vertices(1,1)),(vertices_i(2,1) - vertices(2,1))) * rho_old_time(n) * rho_old_time(i);
end