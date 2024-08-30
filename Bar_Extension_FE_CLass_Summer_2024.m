% Clear the command window
clc

% Clear all the workspace
clear all

% Set the number format to long
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Pre-Processor                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
% (A) Geometric data
%%%%%

% Input the length of the beam
L = 3 ; % <--- input required

% Define the constant corresponding to exponentially decreasing area of
% crossection
% A0 is the area of cross-section at x = 0.
A0 =0.5  ;  % <--- input required

%%%%%
% (B) Material data
%%%%%

% Enter the values of Young's modulus
E =2e11  ; % <--- input required

%%%%%
% (C) Mesh Data
%%%%%
  
% Input the number of elements
N = 320 ; % <--- input required

% Number of nodes per element
m = 2 ; % <--- input required

% Get the total number of nodes. Since we are assuming linear elements so
% the number of nodes are total number of elements + 1
total_nodes = N+1 ;  % <--- input required in terms of N and m

% Degree of freedom per node
dof_per_node = 1 ; % <--- input required 

% Degree of freedom per element
dof_per_element = 2 ; % <--- input required 

% Total number of degrees of freedom are degree of freedom per node times
% the total number of nodes
total_dof = total_nodes*dof_per_node ; % <--- input required in terms of 
                       % dof_per_node and total_nodes

% Construct the connectivity matrix.
connectivity_matrix = zeros(N,dof_per_element) ; % <--- input required in terms of N
for i = 1: N
    connectivity_matrix(i,1) = i     ;  % <--- input required in terms of i
    connectivity_matrix(i,2) = i+1 ;  % <--- input required  in terms of i
end

% Construct the coordinate matrix X
coord_matrix =  linspace(0, L, total_nodes)' ;
 
%%%%%
% (D) Load - both body forces and surface traction 
%%%%%

% The constant corresponding to the axially distributed force
p0 = 0 ; % <--- input required

% Applied boundary load at x = L
P_star = 1e3 ; % <--- input required

%%%%%
% (E) Specify the boundary condition
%%%%%

% This is the node (dof) on which u = 0 and x = 0
iux = 1 ; % <--- input required

% Nodes for which the problem has to be solved
ir = setdiff(1:total_dof,iux)' ; %get clarity

%%%%%
% (F) Define the Gauss quadrature related  parameters
%%%%%

Ngaus = 2 ; % Number of Gass quadrature points in one direction.  
gfa = 1.d0/sqrt(3.d0) ;

gauss_data(1,1) = -gfa ; % First Gauss quadrature point coordinate  
gauss_data (1,2) = 1.d0 ; % First Gauss quadrature point weight
gauss_data (2,1) = gfa  ; % Second Gauss quadrature point coordinate
gauss_data (2,2) = 1.d0 ; % Second Gauss quadrature point weight

 % Total number of Gauss quadrature points.
ngp = Ngaus ;   % <--- input required in terms of Ngaus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Processor                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize the  global stiffness matrix and right side force vector

K_global = zeros(total_dof,total_dof) ; % <--- input required in terms of total_dof 
F_global = zeros(total_dof, 1) ; % <--- input required in terms of total_dof 
x_gauss=zeros(N,ngp);
% Run the loop over all the elements
for e = 1 : N % <--- input required in terms of N 
    
    % Get the connectivity of element e
    con_e = connectivity_matrix(e,:) ; % <--- input required 
    
    % get the coordinate of the element e
    x_local = coord_matrix(con_e,1) ; % <--- input required in terms of con_e
    
    % Calculate the length of the element
    l_e =x_local(2)-x_local(1) ; % <--- input required in terms  x_local
    
    % Initialize the elemental stiffness and right side force vector
    K_elemental = zeros(dof_per_element,dof_per_element) ; % <--- input required in terms of dof_per_element
    F_elemental = zeros(dof_per_element,1) ; % <--- input required in terms of dof_per_element
    
    % Run the Gauss integration loop over ngp Gauss points
    for gp = 1: ngp % <--- input required in terms of ngp
    
        xi = gauss_data(gp,1) ; % xi coordinate of the gp Gauss point
        weight = gauss_data(gp,2) ; % weight of the gp Gauss point
        
        N1 = (1-xi)/2 ; % Shape function at local node 1 % <--- input required
        N2 = (1+xi)/2 ; % Shape function at local node 2 % <--- input required
        
        % Form the elemental shape function matrix
        Shape_function_matrix = [N1 ; N2 ] ; % <--- input required
        
        Be = [-1/2 ; +1/2 ] ;
        
        % Calculate the area at the current Gauss point coordinate
        x =  Shape_function_matrix'*x_local  ;
        x_gauss(e,gp)=x;% <--- input required. Use Shape_function_matrix
        A = A0*exp(-x/L) ; % <--- input required. Assume a suitable function
        p_x = 0 ; % <--- input required. Assume a suitable function
        
        % Compute the local stiffness and right side vector at the current
        % Gauss point 
        
        K_elemental = K_elemental + A*E*Be*Be'*2*weight/l_e ; % <--- input required. 
        F_elemental = F_elemental + p_x*Shape_function_matrix*l_e*0.5*weight; % <--- input required. 
    end
    
    % Assemble the local stiffness and right side vector in the global
    % stiffness matrix and right side vector
    K_global(con_e,con_e) = K_global(con_e,con_e) + K_elemental ; % <--- input required. 
    F_global(con_e,1) = F_global(con_e,1) + F_elemental ; % <--- input required. 
        
end

% Put the force at x = L in the global right side force vector
F_global(total_dof,1) = F_global(total_dof,1) + P_star ; % <--- input required. 

% Apply the boundary conditions on the global stiffness matrix and right
% side force vector
K_hat = K_global(ir,ir) ; % <--- input required. Use ir vector.
F_hat = F_global(ir,1) ; % <--- input required.

% Solve for the unknown displacements
U_ir = K_hat\F_hat ; % <--- input required.

% Form the global displacement matrix
U = zeros(total_dof,1) ; % <--- input required.
U(ir) = U_ir ; % <--- input required.

% Get the displacement of the node at x = L
disp(U(total_dof,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Post-processor                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the displacement
plot(coord_matrix,(U),'-ro') ;
hold on
title(' Plot of displacement of the rod along its length ');
xlabel('     x-axis [m]  -->    ');
ylabel('   Displacement [m]   -->    ');
hold on ;



stresses_gauss = zeros(N, ngp);

% Loop over elements
for e = 1:N
    con_e = connectivity_matrix(e,:) ; % <--- input required

    % get the coordinate of the element e
    x_local = coord_matrix(con_e,1) ; % <--- input required in terms of con_e

    % Calculate the length of the element
    l_e =x_local(2)-x_local(1) ; % <--- input required in terms  x_local
    % Internal force (axial force) in the element
    u1=U(con_e(1));
    u2=U(con_e(2));
   % Loop over Gauss points
    for gp = 1:ngp
        % Gauss point location in the element
        xi = gauss_data(gp,1) ;
        
        % Shape functions and derivative at Gauss point
        N1 = 0.5 * (1 - xi);
        N2 = 0.5 * (1 + xi);
        dN1_dxi = -0.5;
        dN2_dxi = 0.5;
        
        % Jacobian determinant
        J = l_e / 2;  % since dxi/dx = 1 / (2 * J)
        
        % Strain at Gauss point
        epsilon = (dN1_dxi * U(con_e(1)) + dN2_dxi * U(con_e(2)))*(1/(2*J));
        
        % Stress at Gauss point using Hooke's law
        stress_gp = E * epsilon;
        
        % Store stress at Gauss point
        stresses_gauss(e, gp) = stress_gp;
    end
end
 %Plot the Stress
x_1D = reshape(x_gauss.', 1, []);
x_1D=[0,x_1D];
x_1D=[x_1D,3];
stresses_1D=reshape(stresses_gauss.', 1, []);
stresses_1D=[stresses_1D(1),stresses_1D];
stresses_1D=[stresses_1D,stresses_1D(end)];
plot(x_1D,stresses_1D,'-ro') ;
hold on
title(' Plot of stress in the rod along its length ');
xlabel('     x-axis [m]  -->    ');
ylabel('   Stress [N/m2]   -->    ');
hold on ;


disp('Stresses at Gauss points in each element:');
disp(stresses_gauss);
f_stress = fopen('stress_at_GP_mesh_320.txt', 'w');  % Open file for writing
fprintf(f_stress, 'Stresses at Gauss points in each element:\n');
for i=1:N
    fprintf(f_stress, 'Element %d:%f %f \n', i,stresses_gauss(i,:));
end
fclose(f_stress);  % Close the file
f_d = fopen('displacement_along_length_mesh_220.txt', 'w');  % Open file for writing
fprintf(f_d, 'Displacement along the length :\n');
for i=1:total_dof
    fprintf(f_d,'Node %d:%.12e \n',i, U(i));
end
fclose(f_d);


% Task 1
%
% Get the plot of displacement along the length

% Task 2
%
% Get the stresses at the Gauss point and print in a file where rows are
% number of elements and columns are stress values at the Gauss point of
% the elements. Note that since it is a 1D problem there is only one stress
% component.