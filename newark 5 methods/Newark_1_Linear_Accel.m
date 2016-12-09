%% Linear Acceleration (Newark's Method #1)
load('reduced_K_M')
load('C_matrix')

%Declaring variables
endT = 0.13;
T = 0.01;% time
dT = 0.000001;
% final time
Beta = 1/12;
gamma = 1/2;
R1 = zeros(150,1);
R1(149,1) = 100000;
R0 = zeros(150,1);
d = zeros(150,1); dd = zeros(150,1); ddd = M_r\R1;
%ddd = zeros(150,1);
Timestep = 0:dT:endT;
nstep = length(Timestep);
disp = zeros(150,1); velo = disp; acce = velo; time = zeros(150,nstep);
%% For loop required
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

disp(:,i) = real(dn(:));
velo(:,i) = real(ddn(:));
acce(:,i) = real(dddn(:));

d = dn;
dd = ddn;
ddd = dddn;
end

theta = disp(121,:);
dtheta = velo(121,:);
ddtheta = acce(121,:);
hold on
figure(1)
plot(Timestep,theta) 
ylabel('thetaz41 (rad/s)')
xlabel('time(s)')
figure(2)
plot(Timestep,dtheta) 
ylabel('dthetaz41 (rad/s)')
xlabel('time(s)')
figure(3)
plot(Timestep,ddtheta)
ylabel('ddthetaz41 (rad/s)')
xlabel('time(s)')