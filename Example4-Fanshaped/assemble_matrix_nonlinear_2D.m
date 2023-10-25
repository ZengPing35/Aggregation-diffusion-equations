function A=assemble_matrix_nonlinear_2D(P_partition, T_partition, tau, matrix_size, num_control_volume, measure_control_volume, rho_old_n, rho_old_time)
%%用 Up-wind scheme 来组装 $\rho_t = \nabla \cdot \Big[\rho \nabla \big(H'(\rho)+V+W*\rho \big) \Big]$
%% Initialization
A = sparse(matrix_size(1),matrix_size(2)); 
function_discrete_convolution_vector=zeros(num_control_volume,1);
xi_vector=zeros(num_control_volume,1);
grad_H=function_vector_grad_H(rho_old_n);

%% Function convolution; Function grad_H; Function xi
for n = 1:num_control_volume
    vertices=P_partition(:,n);              
    function_discrete_convolution_vector(n)=function_Discrete_convolution(P_partition, num_control_volume, vertices, measure_control_volume, rho_old_n, rho_old_time);
    xi_vector(n) = grad_H(n) + function_V(vertices(1,1), vertices(2,1)) + function_discrete_convolution_vector(n);
end

%% Up-wind scheme
for i = 1:num_control_volume
    A(i, i) = measure_control_volume(i) / tau;  %%rho_t这一项得来的
    xi_i = xi_vector(i);                        
    vertices_A =  P_partition(:,i);             %%当前i控制体的节点坐标   
    [m,n] = find(T_partition == i);             %%n代表控制体i周围三角形的编号
    b = T_partition(:, n);                      %%确定i控制体周围所有三角形节点编号
    c = b(:);               %%矩阵变成向量
    c(c == i) = [];         %%去掉i点
    d = unique(c);          %%去掉重复元素
    num = length(d);  
    for j = 1:num                   %%i控制体周围控制体的个数
        xi_j = xi_vector(d(j));     
        [q,t] = find(b == d(j));    %%判断第j个控制体包含于控制体i周围的哪几个三角形中
        num_vol_in_tri = length(t); 
        vertices_B =  P_partition(:,d(j));  %%第j个控制体包含于第k个三角形的第二个节点坐标
        vertices_AB_mid = (vertices_A + vertices_B)/2;  %%三角形AB边的中点坐标
        edge_AB = sqrt((vertices_B(1) - vertices_A(1))^2 + (vertices_B(2) - vertices_A(2))^2);    %%AB点的距离
        edge_ij = 0;                        %%初始化edge_ij
        for k = 1:num_vol_in_tri
            e = b(:, t(k));                 %%第j个控制体包含于第k个三角形的三角形节点编号
            e(e == i) = [];                 %%第k个三角形除去i节点 
            e(e == d(j)) = [];              %%再去掉第k个三角形的节点
            vertices_C =  P_partition(:,e); %%第j个控制体包含于第k个三角形的第三个节点坐标                   
            vertices_tri_ij = circumcenter(vertices_A(1), vertices_A(2), vertices_B(1), vertices_B(2), vertices_C(1), vertices_C(2));      %%三角形的外心坐标
            
           %% 画出三角形及外心坐标
%             x = [vertices_A(1), vertices_B(1), vertices_C(1), vertices_tri_ij(1)];
%             y = [vertices_A(2), vertices_B(2), vertices_C(2), vertices_tri_ij(2)];
%             plot(x, y, 'r*');
%             axis([-1, 1, -1, 1]);
            
           %% edge:(AB中点到外心线段)--->算的是AB中点到旁边所有三角形外心长度之和           
            edge_ij = edge_ij + sqrt((vertices_tri_ij(1) - vertices_AB_mid(1))^2 + (vertices_tri_ij(2) - vertices_AB_mid(2))^2); %%AB边中点到外心的距离
        end
        upwind = (edge_ij / edge_AB) * (xi_i - xi_j);
        if upwind >= 0
            A(i, i) = A(i,i)+upwind;
        else
            A(i, d(j)) = upwind;
        end
    end
end