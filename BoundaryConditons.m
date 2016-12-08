%%% This Code is used to apply the boundary condition for the Beam Problem
%%% thius obtaining the reduced K and M matrix %%%

load('project3_input_structure','K','M')


%%Applying the Boundary conditions in X-direction
r1 = 7:6:306; 
%%Applying the Boundary Condition in Y-direction
r2 = 8:6:300;
%%Applying the Boundary Condition in Theta_z direction
r3 = 6:6:306;

rows = sort([r1 r2 r3]);
column = rows';
K_r = full(K(rows,column));
M_r = full(M(rows,column));

save('reduced_K_M','K_r','M_r')