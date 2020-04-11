function [stresses,strains,dbar_e] = postprocessing(Xbar,D,GDOF,flag)
% This function determines the approximated stresses and strains from the
% given discrete deflection D
% Inputs:
% Xbar          array holding the absolut position of each node
% D             discrete displacement vector
% GDOF          global DOF 
% flag          structured array containing properties
% Outputs:
% stresses      approximated (discrete) stresses
% strains       approximated (discrete) strains 
% dbar_e        nodal position after deflection D

%% evaluate constitutive matrix C (plane strain)
young = flag.material.emod;
poisson = flag.material.poisson;
C = young/(1.0-poisson^2) ...
 * [ 1.0      poisson  0.0               ;
     poisson  1.0      0.0               ;
     0.0      0.0      (1.0-poisson)/2.0 ];
 
%% calculate nodal position after deflection D
[DOF_e,numele] = size(Xbar);
numnodes = DOF_e/2;
dbar_e = zeros(DOF_e,numele);
for e=1:numele
    for i=1:DOF_e
        idx_i = GDOF(i,e);
        dbar_e(i,e) = D(idx_i);
    end
end
eta = [-1,-1;1,-1;1,1;-1,1];
%% calculate strains and stresses
strains = zeros(3*numnodes,numele);
stresses = zeros(3*numnodes,numele);
for e=1:numele
    xbar = Xbar(:,e);   
    for n=1:(numnodes)
        eta_node = eta(n,:);
        [~,~,B] = evaluate_bilinear_shapefct(xbar,eta_node);
        strain_node = B*dbar_e(:,e);
        strains(n+(n-1)*2:n+(n-1)*2+2,e) = strain_node;
        stress_node = C*strain_node;
        stresses(n+(n-1)*2:n+(n-1)*2+2,e) = stress_node;
        
    end
end
end