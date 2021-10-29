%% EE 514 Inverted Pendulum Setup Script
% K. Clay McKell
close all; clear all; clc;
N = 8.125;
L = 10.75*2.55/100;
M = 28/1e3;
vm0 = 0.69;
g = 9.8;
KcRa = M*g*L/vm0;
Km = 23.133;
DZ = 1.0685;
taum = 0.274;
Jm = M*g*L/vm0/(Km/taum);
J = Jm;
Ka = -4.95;
Kp = -4.84;
KT = -0.1411;
xdat = -14:2:14;
ydat = ones(size(xdat));
%% Testing
A = 0.7774; % Square wave amplitude.
Kf = 0.0273;
Kb = (Jm/taum-Kf)/KcRa;
dz = 0.0367;
% Determined via GA:
fric = [0.978135931332231,1.01874130511885,0.986205707692089,0.922480795286671,0.935772876461489,1.03269473992516,1.01647871181809,0.976391525447778,1.05977676664789,1.01400860537951,1.01388813786157,0.972708707740273,0.985202668303909,0.948810831127879,1.08451371628133];
% S = sim('PlantOL');
% figure
% plot(S.omega)
% disp(tcfinder(S.omega.Time,S.omega.Data));
% Testing
% close all;
% PC = load('PositionCycle.mat');
% plot(PC.PositionCycle.Times+0.1775,PC.PositionCycle.Math1V+2.73793)
% hold on; 
Jp = L^2*M;
J = Jm + Jp/N/N;
% out = sim('PlantOL');
%plot(out.theta.Time,out.theta.Data);
