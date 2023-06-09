
%This is the main script. In the example below it calls BFconfig.m to compute the actual
%parameters used in the simulation.


% run ONE of the following lines
% BFtype=0; Nmic=2; D=0.03; Ts=1/48e3; NT=0; % broadside configuration
BFtype=1; Nmic=3; D=0.015; Ts=1/48e3; NT=0;  % 1st ord end-fire, Dipole
% BFtype=1; Nmic=2; D=0.004; Ts=0; NT=0; % 1st ord end-fire, Dipole
% BFtype=1; Nmic=2; D=0.0214; Ts=1/16e3; NT=1; % 1st ord end-fire, Cardioid 16kHz PCM
% BFtype=1; Nmic=2; D=0.0071; Ts=1/48e3; NT=1; % 1st ord end-fire, Cardioid 48kHz PCM
% BFtype=1; Nmic=2; D=0.004; Ts=1/1.024e6; NT=12; % 1st ord end-fire, Cardioid 1.024MHz PDM
% BFtype=1; Nmic=3; D=0.004; Ts=1/1.024e6; NT=12; % 2nd ord end-fire, Cardioid 1.024MHz PDM
% then run the configuration utility (config=0 makes it interactive)
config=1; % 1 if vars are set here: BFtype(0-1), Nmic(1-8), D[m], Ts[s] and NT
BFconfig; % set position, weight and time delay of each microphone
% run ANY of the desired simulation scripts
BFsim2D; % 2D sim with planar/circular wavefronts at given frequency
BFsim3D; % 3D sim with circular wavefronts at given frequency
BFsimAF; % 2D sim, plot of intensity vs angle of arrival vs frequency
BFsimAFc; % 2D sim, same as BFsimAF but plots a ring instead of a plane
%%%% BFsimA3; % 3D sim, intensity vs angle of arrival at given frequency

% On the opposite in the example below, all parameters are computed manually.
% See figure 18.

% not using BFconfig.m, try one of the following blocks
% parabolic array, capture sound at focal length
Nmic=6; D=0.015/2; Ts=0; NT=0; focalLen=4;
x=([0:Nmic-1]'-(Nmic-1)/2)*D; y=x.^2/(focalLen*D); m=[x,y,zeros(Nmic,1)];
mw=ones(Nmic,1); maxgain=sum(abs(mw)); mt=zeros(Nmic,1);
f=3e3; f1=0.1e3; f2=20e3; config=1; v=343;
BFstr = sprintf('parabolic, %dmics, D=%.1fmm, T=%d*%.2fus',Nmic,D*1000,NT,Ts*1e6);
% linear array, capture sound from specific angle
Nmic=5; D=0.057; Ts=1/48e3; myAngle=pi/4; % any angle from pi/4 to 3*pi/4

x=([0:Nmic-1]'-(Nmic-1)/2)*D; z=zeros(Nmic,1); m=[x,z,z];
mw=ones(Nmic,1); maxgain=sum(abs(mw));
mt=[0:Nmic-1]'*D/v*cos(myAngle); mt=round(mt/Ts)*Ts; % delay is NT*Ts
f=3e3; f1=0.1e3; f2=20e3; config=1; v=343;
BFstr = sprintf('linear, %dmics, D=%.1fmm, T=NT*%.2fus',Nmic,D*1000,Ts*1e6);