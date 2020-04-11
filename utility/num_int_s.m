function [quad_i] = num_int_s(ijk,xbar,flag)
%NUM_INT_S Summary of this function goes here
%   Detailed explanation goes here
% gather Gauss points and Gauss weights
[GP,GW] = Gauss_lookup();
% compute N,B and J of element
if flag.type == "2D-bilinear"
    % N,B evaluated at GP
    i = ijk(1);
    j = ijk(2);
    eta = [GP(i),GP(j)];
    [N,J,B] = evaluate_bilinear_shapefct(xbar,eta);
end
i = ijk(1);
j = ijk(2);

%% evaluate integrant (volume/body load)
b = flag.load.volume;
quad_ib = GW(i)*GW(j)*N'*b*det(J);
%% evaluate integrant (traction)
% necessary components:
[Neta2fixed_,Jeta2fixed_,~] = evaluate_bilinear_shapefct(xbar,[GP(i),-1]);
[Neta1fixed,Jeta1fixed,~] = evaluate_bilinear_shapefct(xbar,[1,GP(j)]);
[Neta2fixed,Jeta2fixed,~] = evaluate_bilinear_shapefct(xbar,[GP(i),+1]);
[Neta1fixed_,Jeta1fixed_,~] = evaluate_bilinear_shapefct(xbar,[-1,GP(j)]);
qeta2fixed_ = sqrt(Jeta2fixed_(1,1)^2+Jeta2fixed_(1,2)^2);
qeta1fixed = sqrt(Jeta1fixed(2,1)^2+Jeta1fixed(2,2)^2);
qeta2fixed = sqrt(Jeta2fixed(1,1)^2+Jeta2fixed(1,2)^2);
qeta1fixed_ = sqrt(Jeta1fixed_(2,1)^2+Jeta1fixed_(2,2)^2);

t = flag.load.traction;

quad_it = GW(i)*Neta2fixed_'*t*qeta2fixed_ +...
          GW(j)*Neta1fixed'*t*qeta1fixed + ...
          GW(i)*Neta2fixed'*t*qeta2fixed + ...
          GW(j)*Neta1fixed_'*t*qeta1fixed_;
%% total:
quad_i = flag.thickness*(quad_ib + quad_it);
end

