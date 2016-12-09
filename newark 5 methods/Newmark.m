<<<<<<< HEAD:newark 5 methods/Newmark.m
%%% Newmark Method
function [displ, velo, acce] = Newmark

load('project3_input_structure.mat','K','M')

%     disp('')
%     disp('Select the Newmark Method to be used');
%     disp('1 - Linear Acceleration Method');
%     disp('2 - Average Acceleration Method');
%     disp('3 - Algorithmically damped Method');
%     disp('4 - Hilber-Hughes-Taylor Method');
%     disp('5 - Fox-Godwin Method');
%     disp('Select Number for the Newmark Method');
%     
%     Method = input('-> ') ;
% 
%     if Method == 1;
%         gamma = 1/2;
%         Beta = 1/6;
%         dT = 0.000001;
%     elseif Method == 2;
%         gamma = 1/2;
%         Beta = 1/4;
%         dT = 0.0001;
%     elseif Method == 3;
%         gamma = 0.55;
%         Beta = (gamma+0.5)^2*0.25;
%         dT = 0.0001;
%     elseif Method == 4;
%         alpha = -0.25;
%         gamma = 0.5*(1-2*alpha);
%         Beta = 0.25*(1-alpha)^2;
%         dT = 0.0001;
%     else Method = 5;
%         gamma = 1/2;
%         Beta = 1/12;
%         dT = 0.000001;
%     end

%%% The K- and M- Matrix generated using the Beam2 code is in 3D. We need
%%% to convert the K and M matrix in 2D and Apply the Boundary Conditions..

[K_r,M_r] = boundary_conditions(K,M);

%%% Calculation of C-Matrix
zeta = 0.02;
[C,fs] = Damping(K_r,M_r,zeta);



%%% Calculating the Time Increment
Beta = 0.25;
gamma = 0.5;
dT = stability(gamma,Beta,zeta,fs);


endT = 0.13;  %%% End Time for Simulation
T = 0.01;     %%% Impulse Time for Force
%2dT = 0.0001;%%% Time Step

%%% Applied Force
R1 = zeros(150,1);
R1(149,1) = 100000;
R0 = zeros(150,1);

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

a = (1/(Beta*dT^2))*M_r+(gamma/(Beta*dT))*C+K_r;
b = R;
c = (M_r*((1/(Beta*dT^2))*d+(1/(Beta*dT))*dd+(1/(2*Beta)-1)*ddd));
dl = (C*((gamma/(Beta*dT))*d+(gamma/Beta-1)*dd+(gamma/Beta-2)*(dT/2)*ddd));
dn = a\(b+c+dl);


ddn= (gamma/(Beta*dT))*(dn-d)- ((gamma/Beta)-1)*dd - dT*((gamma/(2*Beta))-1)*ddd;

dddn = ((1/(Beta*dT^2))*(dn- d-dT*dd)...
    -((1/(2*Beta))-1)*ddd);

displ(:,i) = real(dn(:));
velo(:,i) = real(ddn(:));
acce(:,i) = real(dddn(:));

d = dn;
dd = ddn;
ddd = dddn;
end

theta = displ(121,:);
dtheta = velo(121,:);
ddtheta = acce(121,:);
hold on;
figure()
plot(Timestep,theta)
title('Displacement')
ylabel('\theta_{z41} (rad/s)')
xlabel('time(s)')
figure()
plot(Timestep,dtheta)
title('Velocity')
ylabel('d\theta_{z41} (rad/s)')
xlabel('time(s)')
figure()
plot(Timestep,ddtheta)
title('Acceleration')
ylabel('dd\theta_{z41} (rad/s)')
xlabel('time(s)')