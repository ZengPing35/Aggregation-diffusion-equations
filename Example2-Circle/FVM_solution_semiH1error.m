function result=FVM_solution_semiH1error(rho_exact, rho, P_partition, T_partition, num_control_volume)
%Numerically compute a norm error of FE solution on the whole Circle domain.
%rho: the values of the FE solution at all nodes of FE in the whole domain. These values are in 1D index of nodes of FE.
%accurate_function: the accurate function in the error.

%% Parameters
result1=0;%numerator of a fraction
result2=0;%the denominator of a fraction

%% Go through all elements and accumulate the error on them.
for i = 1:num_control_volume
    vertices_A =  P_partition(:,i);             %%当前i控制体的节点坐标 
    rho_i = rho(i);
    rho_exact_i = rho_exact(i);
    [m,n] = find(T_partition == i);             %%n代表控制体i周围三角形的编号
    b = T_partition(:, n);                      %%确定i控制体周围所有三角形节点编号
    c = b(:);               %%矩阵变成向量
    c(c == i) = [];         %%去掉i点
    d = unique(c);          %%去掉重复元素
    num = length(d);  
    for j = 1:num                   %%i控制体周围控制体的个数
        [q,t] = find(b == d(j));    %%判断第j个控制体包含于控制体i周围的哪几个三角形中
        num_vol_in_tri = length(t); 
        vertices_B =  P_partition(:,d(j));  %%第j个控制体包含于第k个三角形的第二个节点坐标
        rho_j = rho(d(j));
        rho_exact_j = rho_exact(d(j));
        vertices_AB_mid = (vertices_A + vertices_B)/2;  %%三角形AB边的中点坐标
        edge_AB = sqrt((vertices_B(1) - vertices_A(1))^2 + (vertices_B(2) - vertices_A(2))^2);    %%AB点的距离
        edge_ij = 0;                        %%初始化edge_ij
        for k = 1:num_vol_in_tri
            e = b(:, t(k));                 %%第j个控制体包含于第k个三角形的三角形节点编号
            e(e == i) = [];                 %%第k个三角形除去i节点 
            e(e == d(j)) = [];              %%再去掉第k个三角形的节点
            vertices_C =  P_partition(:,e); %%第j个控制体包含于第k个三角形的第三个节点坐标                   
            vertices_tri_ij = circumcenter(vertices_A(1), vertices_A(2), vertices_B(1), vertices_B(2), vertices_C(1), vertices_C(2));      %%三角形的外心坐标
            
           %% edge:(AB中点到外心线段)--->算的是AB中点到旁边所有三角形外心长度之和           
            edge_ij = edge_ij + sqrt((vertices_tri_ij(1) - vertices_AB_mid(1))^2 + (vertices_tri_ij(2) - vertices_AB_mid(2))^2); %%AB边中点到外心的距离
        end
        rho_exact_ij = rho_exact_i - rho_exact_j;
        rho_numerical_ij = rho_i - rho_j;
        result1 = result1 + (abs(rho_exact_ij - rho_numerical_ij))^2 * edge_ij / edge_AB;
        result2 = result2 + (abs(rho_exact_ij))^2 * edge_ij / edge_AB;
    end
end

%% semi-norm
result1=sqrt(result1);
result2=sqrt(result2);
result=result1/result2;