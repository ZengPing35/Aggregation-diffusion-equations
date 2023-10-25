function r=function_Discrete_convolution(P_partition, number_of_control_volume, vertices_i, measure_control_volume, rho_old_n, rho_old_time)

r=0;
for n = 1:number_of_control_volume
    vertices = P_partition(:,n);  
    r = r + (1/2) * measure_control_volume(n) * function_W((vertices_i(1,1) - vertices(1,1)),(vertices_i(2,1) - vertices(2,1))) * (rho_old_n(n) + rho_old_time(n));
end