function [quad_i] = num_int_k(ijk,xbar,flag)
% This function performs Gauss quadrature to compute k_e of element e.
% Inputs:
% ijk array holding the indices for the Gauss quadrature
% xbar  array holding the absolut position of each node of element e
% flag  structured array containing properties
% Outputs:
% quad_i array of size [e_DOF,e_DOF]

% gather Gauss points and Gauss weights
[GP,GW] = Gauss_lookup();
% evaluate constitutive matrix C (plane strain)
young = flag.material.emod;
poisson = flag.material.poisson;
C = young/(1.0-poisson^2) ...
 * [ 1.0      poisson  0.0               ;
     poisson  1.0      0.0               ;
     0.0      0.0      (1.0-poisson)/2.0 ];
 
% compute B,J of element e
if flag.type == "2D-bilinear"
  % N,B evaluated at GP
    i = ijk(1);
    j = ijk(2);
    eta = [GP(i),GP(j)];
    [~,J,B] = evaluate_bilinear_shapefct(xbar,eta);
end

% evaluate integrant
i = ijk(1);
j = ijk(2);
quad_i = GW(i)*GW(j)*B'*C*B*det(J);
end

