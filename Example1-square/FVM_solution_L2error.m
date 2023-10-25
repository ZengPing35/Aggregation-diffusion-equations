function result=FVM_solution_L2error(num_control_volume, measure_control_volume, rho_exact, rho)
%Numerically compute a norm error of FE solution on the whole Circle domain.
%rho_exact: the accurate function in the error.
%rho: numerical solution.
%Relative error

result1=0;  %%numerator of a fraction
result2=0;  %%the denominator of a fraction

%% Go through all elements and accumulate the error on them.
for n = 1:num_control_volume
    result1 = result1 + abs(rho_exact(n)-rho(n))^2 * measure_control_volume(n);
    result2 = result2 + abs(rho_exact(n))^2 * measure_control_volume(n);
end
 
%% L2-norm
result1=sqrt(result1);
result2=sqrt(result2);
result=result1/result2;