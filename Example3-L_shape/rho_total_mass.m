function r=rho_total_mass(rho_h, measure_control_volume, rho_num_control_volume)

%% Initialization
r = 0;

%%
for i = 1:rho_num_control_volume
    r=r+rho_h(i)*measure_control_volume(i);   
end
