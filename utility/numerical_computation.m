function [K_e,S0_e] = numerical_computation(Xbar,flag)
% This function calculates the element stiffness matrix and the element load vector for each element.
% Inputs:
% Xbar  array holding the absolut position of each node
% flag  structured array containing properties
% Outputs:
% K_e  array holding the element stiffness matrix of each element
% S0_e  array holding the element load vector of each element
K_e = {};
S0_e = {};
%% numerical integration
if flag.type == "2D-bilinear"
    n = 2;
    numGP = ceil((n+1)/2);
    for e = 1:flag.numele
        xbar = Xbar(:,e);
        k_e = zeros(8,8);
        s0_e = zeros(8,1);
        for j = 1:numGP
            for i = 1:numGP
                ijk = [i,j];
                % element stiffness matrix: integrant: B'CB|J| 
                k_e = k_e + num_int_k(ijk,xbar,flag);
                % element load vector: 
                s0_e = s0_e + num_int_s(ijk,xbar,flag);
            end
        end
        K_e{e} = k_e;
        S0_e{e} = -s0_e;    % Don't forget to negate!
    end
K_e = cell2mat(K_e);
S0_e = cell2mat(S0_e);
end