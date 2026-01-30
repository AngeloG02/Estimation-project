%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANT-X SIMULATOR - MAIN                                                  %
% Authors:  Mattia Giurato (mattia.giurato@polimi.it)                     %
%           Paolo Gattazzo (paolo.gattazzo@polimi.it)                     %
% Date: 13/12/2017                                                        %
% Adapted to ANT-X 2DoF by:  Salvatore Meraglia (salvatore.meraglia@polimi.it)%
% Date: 22/12/2022                                                        %
%
% Further modified to include structure three-state identified longitudinal model
% 06/01/23 ML
% 
% Further modified to pick outputs with measurement error
% 03/01/24 ML
% 06/01/25 ML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearvars;

% addpath('datasets','common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');

%% Model parameters

% Initial model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)

Xu=-0.1068;


Xq=0.1192;

Mu=-5.9755;

Mq=-2.6478;

Xd=-10.1647;

Md=450.71;

A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];

B=[Xd; Md; 0];

C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 

D=[0; 0; 0; Xd];

% Noise

%noise.Enabler = 0;
noise.Enabler = 1;

noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]

noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]

noise.acc_stand_dev = noise.Enabler * 0.5;                               %[m/s^2]

noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);                   %[rad/s]

% Delays

delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

%% Load controller parameters

parameters_controller                    

%% M injection example (sweeep: first column time vector, secondo column time history of pitching moment) 

load ExcitationM

SetPoint=[0,0];

%% Values selected

t=ExcitationM(:,1);

simulation_time=t(end)-t(1);

%% Launch SIMULATOR

out = sim("Simulator_Single_Axis");


%% Delete temporary files

if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end

%% END OF CODE