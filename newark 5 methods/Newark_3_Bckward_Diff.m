%% Backward Difference (Newark's Method #3)

%Declaring variables
T = 0;              % time
dT = 0.01;          % final time
beta = 1;
gamma = 3/2;
R = 100; % Force arbittary value taken for now
[M] = [];           % Mass Matrix
[K] = [];           % Stiffness matrix
[C] = [];           % Damping Matrix
[D] = [];           % Displacement Matrix
[V] = [];           % Velocity Matrix
[A] = [];           % Acceleration Matrix

%% For loop required

Disp_Future = (inv((1/(beta*dT^2))*[M] + [K]))*(R*(T+dT) + [M]*((1/(beta*dT^2))*[D] + (1/(beta*dT))*[V] + ((1/(2*beta))-1)*[A])
    + [C]*((gamma/(beta*dT))*[D] + ((gamma/beta)-1)*[V] + ((gamma/beta)-2)*(DT/2)*[A]))

Velo_Future = ((gamma/(beta*dT))*(Disp_Future- [D])) - (gamma/beta)-1)*[V] - dT*((gamma/(2*beta))-1)*[A]

Accl_Future = ((1/(beta*dT^2))*(Disp_Future- [D]-dT*[V])-((1/(2*beta))-1)*[A])