function [K,F] = assembly(k_e_all,s0_e_all,Xbar,GDOF)
% This function assembles the elementwise quantities k_e/s0_e to the global matrix/vector K/F.
% Inputs:
% k_e_all       array of elementwise stiffness matrices
% s0_e          array of elementwise load vectors 
% Xbar          array holding the absolut position of each node
% GDOF          global DOF
% Outputs:
% K             global stiffness matrix of our (discretized) structure
% F             global load vector of our (discretized) structure
%% Some pre-processing
[DOF_e,numele] = size(Xbar);
DOF = max(max(GDOF));
EDOF = [1:DOF]';
EDOF = repmat(EDOF,1,numele);
K = zeros(DOF,DOF);
F = zeros(DOF,1);
%% Assembly of K
for i=1:DOF
    for j = 1:DOF
        Kij = 0;
        for e=1:numele
            k_e = k_e_all(:, 1+DOF_e*(e-1):DOF_e*(e-1)+DOF_e);
            % does element contribute to global DOF ij?
            if ismember(i,GDOF(:,e)) && ismember(j,GDOF(:,e))
                % -> yes: add corresponding entry of k_e to global K
                idx_i = GDOF(:,e) == i;
                idx_j = GDOF(:,e) == j;
                idx_i_loc = EDOF(idx_i,e);
                idx_j_loc = EDOF(idx_j,e);
                Kij = Kij + k_e(idx_i_loc,idx_j_loc);
            end
        end
        K(i,j) = Kij;
    end
end
%% Assembly of F
for i = 1:DOF
    Fi = 0;
    for e=1:numele
        s0_e = s0_e_all(:,e);
        % does element contribute to global DOF i?
        if ismember(i,GDOF(:,e))
            % -> yes: add corresponding entry of s0_e to global F
            idx_i = GDOF(:,e) == i;
            idx_i_loc = EDOF(idx_i,e);
            Fi = Fi - s0_e(idx_i_loc);
        end
    end
    F(i) = Fi;
end
end