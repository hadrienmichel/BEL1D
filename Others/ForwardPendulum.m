function [Y] = ForwardPendulum(times,param)%l,h,M)
% FORWARDPENDULUM is a fucntion that computes the forward model for the
% simple pendulum experiment. It takes as arguments:
%       - times: the times of measurements
%       - param = [l, h, M]
%           - l: the length of the pendulum (in meters)
%           - h: the initial heigth of the mass (in meters)
%           - M: the mass of the pendulum (in kg)
%       
% The function returns a vector with the Y position of the pendulum at the 
% times in the vector 'times'.

%tic;
l = param(1);
h = param(2);
M = param(3);
%% 1) Tranforming the h and l values into the angular position of the
%     pendulum at t_0

H = 10; %Height of the pivot point, fixed by the experimental design
tmp = H - h;% Removing the component above the ground
theta_0 = acos(tmp/l);% Angular position of the pendulum at t_0
omega_0 = sqrt(9.81/l);% Value of the constant omega_0

%% 2) Definition of the position with time:
f = @(t,in2,param1)[in2(2,:);-param1.^2.*sin(in2(1,:))];
% The equation is solving the simple pendulum with
%       t = time
%       in2 = [angular position; angular cvelocity]
%       param1 = omega_0;
% The equation returns x:
%       x(:,1) = angular position
%       x(:,2) = angular velocity

x0 = [theta_0; 0];
%t0 = 2*pi()/omega_0;
tInit = 0;
if(size(times,2) ~= 1)
    tMeasure = times;
else
    tMeasure = times';
end
[t, x] = ode45(@(t,in2) f(t,in2,omega_0), [tInit tMeasure], x0);
% t is the time vector and x is the vector containing the values of 
%       x(:,1) angular position
%       x(:,2) velocity

%% 3) Finding the values that are the closests to X_observed
Y = ones(length(times),1).*10 - cos(x(2:end,1)).*l;
%toc

end
