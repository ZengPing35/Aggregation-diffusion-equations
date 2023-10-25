function result=FE_solution_error_triangle_index(rho, accurate_function, t, P_partition, T_partition, basis_type, basis_index, derivative_degree_x, derivative_degree_y, Gauss_point_number)
%Numerically compute a norm error of FE solution on the whole circle domain.
%rho: the values of the FE solution at all nodes of FE in the whole domain. These values are in 1D index of nodes of FE.
%accurate_function: the accurate function in the error.
%basis_type: the type of the FE.
%basis_type=0:2D constant FE.
%basis_type=1:2D Lagrange linear FE.
%basis_type=2:2D Lagrange quadratic FE.
%basis_type=10:2D Crouzeix-Raviart FE.
%derivative_degree_x:the derivative degree of the FE solution with respect to x.
%derivative_degree_y:the derivative degree of the FE solution with respect to y.
%Gauss_point_number: the number of the Gauss points of the Gauss quadrature we want to use.
%vertices: the coordinates of the vertices of a triangular element.
%rho_h_local: the values of the FE solution at the nodes of FE in a triangular element.

number_of_elements = size(T_partition,2);            %%三角形个数
T_basis=T_partition;
[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

result1=0;
result2=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    vertices=P_partition(:,T_partition(:,n));
    rho_local=rho(T_basis(:,n));
    result1=result1+Gauss_quad_2D_error_triangle_time(rho_local,accurate_function,t,vertices,basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,derivative_degree_x,derivative_degree_y,1);
    result2=result2+Gauss_quad_2D_error_triangle_time(rho_local,accurate_function,t,vertices,basis_index,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,derivative_degree_x,derivative_degree_y,0);
end
result1=sqrt(result1);
result2=sqrt(result2);
result=result1/result2;