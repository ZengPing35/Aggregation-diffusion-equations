function measure_control_volume = function_measure_control_volume(P_partition, T_partition, num_control_volume)

%% Initialization
measure_control_volume = zeros(num_control_volume, 1);

%% 
for i = 1:num_control_volume
    [m,n] = find(T_partition(1:3, :) == i);    
    num_triangle = length(n);
    vertices_A =  P_partition(:,i); %%当前i控制体的节点坐标
    measure_control_volume_i = 0;   %%控制体i的面积
    
    for j = 1:num_triangle
        a = T_partition(1:3, n(j));     %%确定i控制体周围的第j个三角形节点编号
        a(a == i) = [];                 %%得到i控制体周围的第j个三角形的另外两个节点编号        
        vertices_B =  P_partition(:,a(1)); %%当前i控制体周围的第j个三角形的第一个节点坐标
        vertices_C =  P_partition(:,a(2)); %%当前i控制体周围的第j个三角形的第二个节点坐标
        vertices_AB_mid = (vertices_A + vertices_B)/2;  %%三角形AB边的中点坐标
        vertices_AC_mid = (vertices_A + vertices_C)/2;  %%三角形AC边的中点坐标
        vertices_tri_ij = circumcenter(vertices_A(1), vertices_A(2), vertices_B(1), vertices_B(2), vertices_C(1), vertices_C(2));      %%三角形的外心坐标
        %% 画出三角形及外心坐标
%             x = [vertices_A(1), vertices_B(1), vertices_C(1), vertices_tri_ij(1)];
%             y = [vertices_A(2), vertices_B(2), vertices_C(2), vertices_tri_ij(2)];
%             plot(x, y, 'r*');
%             axis([-1, 1, -1, 1]);
        %%
        edge_A_ABmid = sqrt((vertices_AB_mid(1) - vertices_A(1))^2 + (vertices_AB_mid(2) - vertices_A(2))^2);    %%AB边中点到i控制体的距离(这里只计算的AB边长度的一半)
        edge_A_ACmid = sqrt((vertices_AC_mid(1) - vertices_A(1))^2 + (vertices_AC_mid(2) - vertices_A(2))^2);    %%AC边中点到i控制体的距离(这里只计算的AC边长度的一半)
        edge_midAB_circ = sqrt((vertices_AB_mid(1) - vertices_tri_ij(1))^2 + (vertices_AB_mid(2) - vertices_tri_ij(2))^2);    %%AB边中点到外心的距离
        edge_midAC_circ = sqrt((vertices_AC_mid(1) - vertices_tri_ij(1))^2 + (vertices_AC_mid(2) - vertices_tri_ij(2))^2);    %%AC边中点到外心的距离
        measure_control_volume_i =  measure_control_volume_i + 1/2 * edge_A_ABmid * edge_midAB_circ + 1/2 * edge_A_ACmid * edge_midAC_circ;   %%控制体i的测度
    end
    measure_control_volume(i) = measure_control_volume_i;      %%i控制体的测度    
end
