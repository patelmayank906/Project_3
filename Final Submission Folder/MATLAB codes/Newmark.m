%%% Newmark Method
%%% Will only work for the problem discussed in the textbook. This is not a
%%% generalized code


%%% d = displacement at current time step
%%% dn = displacement at next time step
%%% dd = velocity at current time step
%%% ddn = Velocity at next time step 
%%% ddd = acceleration at current time step
%%% dddn = acceleration at next time step


function [displ, velo, acce] = Newmark

load('project3_input_structure.mat','K','M');

    disp('');
    disp('Select the Newmark Method to be used');
    disp('1 - Linear Acceleration Method');
    disp('2 - Average Acceleration Method');
    disp('3 - Algorithmically damped Method');
    disp('4 - Hilber-Hughes-Taylor Method');
    disp('5 - Fox-Godwin Method');
    disp('Select Number for the Newmark Method');
    
    Method = input('-> ') ;

    if Method == 1;
        gamma = 1/2;
        Beta = 1/6;
        dT = 0.000001;
    elseif Method == 2;
        gamma = 1/2;
        Beta = 1/4;
        dT = 0.0001;
    elseif Method == 3;
        gamma = 0.55;
        Beta = (gamma+0.5)^2*0.25;
        dT = 0.0001;
    elseif Method == 4;
        alpha = -0.25;
        gamma = 0.5*(1-2*alpha);
        Beta = 0.25*(1-alpha)^2;
        dT = 0.0001;
    else
        gamma = 1/2;
        Beta = 1/12;
        dT = 0.000001;
    end

%%% The K- and M- Matrix generated using the Beam2 code is in 3D. We need
%%% to convert the K and M matrix in 2D and Apply the Boundary Conditions..

[K_r,M_r] = boundary_conditions(K,M);

%%% Calculation of C-Matrix
zeta = 0.02;
[C,fs] = Damping(K_r,M_r,zeta);

%%% Calculating the Time Increment
% Beta = 0.25;
% gamma = 0.5;
% dT = stability(gamma,Beta,zeta,fs);


endT = 0.13;  %%% End Time for Simulation
T = 0.01;     %%% Impulse Time for Force
%dT = 0.0001;%%% Time Step

%%% Applied Force
R1 = zeros(150,1);
R1(149,1) = 100000; %%% Force applied upto Time T
R0 = zeros(150,1); %%% Force applied after Time T

%%% Initial Conditions
d = zeros(150,1); dd = zeros(150,1); ddd = M_r\R1;

%%% TimeStep Calculation
Timestep = 0:dT:endT;
nstep = length(Timestep);

%%% Preallocation of Matrices
displ = zeros(150,1); velo = displ; acce = velo;


%%% Loop to Calculate Displacement, Velocity and Acceleration

for i = 1:nstep;
    step = Timestep(i);
if Timestep(i) <= T
    R = R1;
else
    R = R0;
end

%%% Calculate the displacement for the Next Time Step
a = (1/(Beta*dT^2))*M_r+(gamma/(Beta*dT))*C+K_r;
b = R;
c = (M_r*((1/(Beta*dT^2))*d+(1/(Beta*dT))*dd+(1/(2*Beta)-1)*ddd));
dl = (C*((gamma/(Beta*dT))*d+(gamma/Beta-1)*dd+(gamma/Beta-2)*(dT/2)*ddd));
dn = a\(b+c+dl);

%%% Calculate the Velocity for the Next Time Step
ddn= (gamma/(Beta*dT))*(dn-d)- ((gamma/Beta)-1)*dd - dT...
    *((gamma/(2*Beta))-1)*ddd;

%%% Calculate the Acceleration for the Next Time Step
dddn = ((1/(Beta*dT^2))*(dn- d-dT*dd)...
    -((1/(2*Beta))-1)*ddd);


%%% Saving the Value at Iteration in a Matrix
displ(:,i) = real(dn(:));
velo(:,i) = real(ddn(:));
acce(:,i) = real(dddn(:));

%%% Reset the Value for the Next Iteration
d = dn;
dd = ddn;
ddd = dddn;
end
save('model_parameters','velo','displ','acce'); % Store desplacement
                                                % velcity and acceleration 
                                                % value in a .mat file
theta = displ(121,:);
dtheta = velo(121,:);
ddtheta = acce(121,:);

%%% Plot Figures
hold on; 
figure();grid on
plot(Timestep,theta)
title('Displacement')
ylabel('\theta_{z41} (rad/s)')
xlabel('time(s)')

figure();grid on
plot(Timestep,dtheta)
title('Velocity')
ylabel('d\theta_{z41} (rad/s)')
xlabel('time(s)')

figure();grid on
plot(Timestep,ddtheta)
title('Acceleration')
ylabel('dd\theta_{z41} (rad/s)')
xlabel('time(s)')