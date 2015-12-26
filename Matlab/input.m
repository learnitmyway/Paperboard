% Input for simplified example

clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Discrete (pseudo) time intervals
t_data = [0 5 10 15 20 25 30 35]; % [0 5 10 15 20 25 30 35]

% Strain component 
eps33_data = [0 5e-3 0 -5e-3 0 5e-3 0 -5e-3 0 5e-3]; % [0 5e-3 0 -5e-3 0 5e-3 0 -5e-3 0 5e-3]
eps13_data = [0 0 0 0 0 0 0 0]; % 
eps23_data = [0 0 0 0 0 0 0 0]; % [0 0 0 0 0 0 0 0]

% Number of increments per time interval
N = 50; 

% Linear interpolation
[time,eps33] = linearInterpolation(t_data,eps33_data,N);
[time,eps13] = linearInterpolation(t_data,eps13_data,N);
[time,eps23] = linearInterpolation(t_data,eps23_data,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation (t = 0)
stran = zeros(6,1);
dstran = zeros(6,1);
statev = zeros(13,1);
stress = zeros(6,1);
ddsdde = zeros(6,6);

% Calculation

for n = 1 : (length(time) - 1)
    
    % Strain (3x1)
    stran(3,n+1) = eps33(n+1); 
    stran(5,n+1) = 2*eps13(n+1);
    stran(6,n+1) = 2*eps23(n+1);
    
    % Strain increment
    dstran(:,n+1) = stran(:,n+1)-stran(:,n);
    
    % Initialise array within loop
    stress(:,n+1) = stress(:,n);
    statev(:,n+1) = statev(:,n);
    ddsdde(:,:,n+1) = ddsdde(:,:,n);
    
    % umat
    [stress(:,n+1),statev(:,n+1),ddsdde(:,:,n+1)] = umat(stress(:,n+1),dstran(:,n+1),stran(:,n+1),statev(:,n+1),ddsdde(:,:,n+1),time(n+1));

end

% Analytical solution

eps33_an = [1/1000, 1/200, 1/440, -1/200]; % 1/1000, 1/200, 1/440, -1/200
S33_an = [1, 15/11, -15/11, -245/121]; % 1, 15/11, -15/11, -245/121
gam23_an = []; % 1/500, 1/100, 1/300, -1/100
S23_an = []; %  1, 5/3, -5/3, -25/9

% plot
plotoutput
