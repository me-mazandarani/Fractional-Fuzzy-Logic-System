%% The reference
% Mehran Mazandarani, Li Xiu, Fractional Fuzzy Inference System:
% The New Generation of Fuzzy Inference Systems,
% IEEE Access, Vol. 8, pp. 126066-126082, 2020. 
% DOI: 10.1109/ACCESS.2020.3008064 
%  Link: https://ieeexplore.ieee.org/document/9136681

% If there is any question contact the authors of the above paper. 
% Date of issue 2020-12-10
% Before using the this file and FFIS.m read carefully the help file,
% help_FFIS.pdf
%% The example of control of inverted pendulum system using FFIS
% This program is the first version of the function FFIS.m (FFIS_201210)
% The algorith of this function has not been written in an optimal manner.
% Thus, naturally, it may take more time than that one may expect for 
% getting the output.
% Based on the settings in this example, it takes almost 7 minutes to get 
% the output.

%% Initialize the Fuzzy system structure (Initialization Section) 
clc
clear
% close all

fis=readfis('fis');

%% Fractional Indices
%     Fids={value_1,form_1, value_2,form_2,...,value_m,form_m}
% Read the help_FFIS.pdf for more information about the cell arraye Fids. 

% fis.Outputs.MembershipFunctions

% In this example the following indices have been considered arbitrary. 
 Fids={0.5, 'a', 0.5, 'a', 1,'a', 0.5, 'b', 0.5,'b'};

% For Mamdani's FIS
%  Fids={1, 'a', 1, 'b', 1,'a', 1, 'b', 1,'b'}; 

%% Model (The main section)




T=3;
n=2000;
tspan=linspace(0,T,n+1);
h=tspan(2)-tspan(1);

% % This is the gain of control signal based on the control structure.
Kgain=220;

% % % % Pendulum parameters
g=9.8;
m=2;
M=8;
l=2;
% x1 is the angle of the pendulum and is used as the error
% x2 is the derivative of the x1, i.e. the derivative of the error
% % % % % % % % % % % % % % % % % %
% x1=theta;
% x2=theta_dot;
x1=zeros(1,n+1);
x2=zeros(1,n+1);
u=zeros(1,n);
x1(1)=0.3;  % the initial condition
x2(1)=0.1;  % the initial condition
a=1/(M+m);

% Model part
for cnt=1:n

Inputs=[x1(cnt); x2(cnt)];

% % % the control signal 
u(cnt)=FFIS(fis,Inputs,Fids);

u(cnt)=Kgain*u(cnt);

% %  The model has been considered in a descrete form by the use of 
% %  forward approximation of the derivative definition.

x1(cnt+1)=x1(cnt)+h*x2(cnt);
k1=g*sin(x1(cnt))-a*m*l*x2(cnt)^2*sin(2*x1(cnt))/2-a*cos(x1(cnt))*u(cnt);
k2=4*l/3-a*m*l*cos(x1(cnt))^2;
x2(cnt+1)=x2(cnt)+h*k1/k2;

end


figure(1)
plot(tspan(1:numel(x1)),x1,'LineWidth',2)
grid on
hold on
figure(2)
plot(tspan(1:numel(u)),u/Kgain,'LineWidth',2)
grid on
hold on






