function [K,F] = enforceDBC(K,F,DBC)
% This function introduces boundary conditions to the system.
% Inputs:
% K             global stiffness matrix of our (discretized) structure
% F             global load vector of our (discretized) structure
% DBC           array that consists of all fixed DOF due to DBC          
% Outputs:
% K             global stiffness matrix of our (discretized) structure
% F             global load vector of our (discretized) structure
K(DBC,:) = [];
K(:,DBC) = [];
F(DBC) = [];
end

