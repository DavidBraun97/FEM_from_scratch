function [D] = solveFEM(K_red,F_red,DBC)
% This function introduces boundary conditions to the system.
% Inputs:
% K_red         reduced global stiffness matrix of our (discretized) structure
% F_red         reduced global load vector of our (discretized) structure
% DBC           array that consists of all fixed DOF due to DBC          
% Outputs:
% D             discrete displacement vector 

%% compute solution to reduced problem
D_red = K_red\F_red;
%% Reassemble to solution D (including the fixed/prescribed deflections)
n = length(F_red) + length(DBC);
D = zeros(n,1);
j = 1;
for i=1:n
    if ~ismember(i,DBC)
        D(i) = D_red(j);
        j = j+1;
    end   
end
end

