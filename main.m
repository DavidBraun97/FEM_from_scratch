% Finite Element Method (FEM) from scratch
%
% 05/20 David Braun
% TU Muenchen
%
clear,clc,%close all
disp('TO DO: traction')

% load corresp. files
addpath('utility');
%% Geometry
% provide the rectangular geometry:
% b//a = 0;
%    C| ----a----
%    L| |      |
%    A| |      |
%    M| b      b
%    P| |      |
%    E| |      |
%    D| ----a----
a = 1; 
b = 3;
visualizegeometry(a,b)
% loads:
Fg_vol = (9.81);
load = struct('volume',[0; -Fg_vol],...
              'traction',[0;0]);
%% Discretization
% relative size of element compared to smaller edge length 
h_rel = 1/3;
T = 1;

mat = struct('name', 'Hooke', ...
             'emod', 10.e6, ...
             'poisson', 0.0);
% mat = struct('name', 'Hooke', ...
%              'emod', 1, ...
%              'poisson', 0.0);
flag = struct('type', '2D-bilinear',  'thickness', T,  'material', mat,'numele',[],'load',load);
%   C|  ----a----
%   L|  |. ..    |
%   A|  |.       |
%   M|  b. ..    b
%   P|  |4-3  8-7|
%   E|  |1-2  5-6|
%   D|  ----a----
% Xbar contains abs. postiton of x1, x2 ...
[Xbar,EDOF,GDOF,flag] = mashing(h_rel,a,b,flag);
%% elementwise computation of stiffness matrix k^(e) and load vector s^0(e)
[k_e,s0_e] = numerical_computation(Xbar,flag);
%% Assembly to global stiffness matrix K and global force vector F
[K,F] = assembly(k_e,s0_e,Xbar,GDOF);
%% Apply Dirichlet Boundary Conditions & external point loads
% locked/supported global DOF
DBC = [1,2,3,4, 9,10,11,12, 17,18]';    % DBC consists of all fixed DOF
ext_force = [77,-10000;...
             79,-10000;...
             75,-10000];
F(ext_force(:,1)) = F(ext_force(:,1)) + ext_force(:,2);

[K_red,F_red] = enforceDBC(K,F,DBC);
%% solve LGS KD=F for D
D = solveFEM(K_red,F_red,DBC);
%% post-processing -> compute resulting stress & strain due to deflection u 
[stresses,strains,dbar_e] = postprocessing(Xbar,D,GDOF,flag);
%% visualize results
visualizeresults(Xbar,dbar_e,GDOF,stresses,strains)