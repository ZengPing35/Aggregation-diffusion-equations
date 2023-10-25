function r = function_free_energy(P_partition, num_control_volume, rho_old_time, measure_control_volume)

%% Initialization
r = 0;
function_discrete_convolution_vector=zeros(num_control_volume,1);
H=function_vector_H(rho_old_time);
xi_vector=zeros(num_control_volume,1);

%% Function convolution; Function grad_H; Function xi
for n = 1:num_control_volume
    vertices=P_partition(:,n);              
    function_discrete_convolution_vector(n) = function_Discrete_convolution_free_energy(P_partition, vertices, num_control_volume, measure_control_volume, rho_old_time, n);
    xi_vector(n) = H(n) + (rho_old_time(n) * function_V(vertices(1,1), vertices(2,1))) + function_discrete_convolution_vector(n);
end


%% Free energy
for i = 1:num_control_volume    
  r = r + measure_control_volume(i) * xi_vector(i);
end
