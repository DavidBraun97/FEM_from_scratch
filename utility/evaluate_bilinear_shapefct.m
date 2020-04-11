function [N,J,B] = evaluate_bilinear_shapefct(xbar,eta)
% This function provides N,J and B of element e evaluated at eta.
% Inputs:
% xbar  array holding the absolut position of each node of element e
% eta   array [eta1,eta2] coordinates of interest
% Outputs:
% N     matrix of shape functions of element e evaluated at eta
% J     Jacobian matrix of element e evaluated at eta
% B     B-operator of element e evaluated at eta

%% evaluate shape functions
eta1 = eta(1);  % coordinate eta_1
eta2 = eta(2);  % coordinate eta_2

N1 = 1./4.*(1.0-eta1)*(1.0-eta2);  % N^1
N2 = 1./4.*(1.0+eta1)*(1.0-eta2);
N3 = 1./4.*(1.0+eta1)*(1.0+eta2);
N4 = 1./4.*(1.0-eta1)*(1.0+eta2);
%% compute N matrix
N = [N1,0,N2,0,N3,0,N4,0;...
     0,N1,0,N2,0,N3,0,N4];

% evaluate derivative shape functions
N1deta1 = 1./4.*(-1.0)*(1.0-eta2);  % N^1_{,eta1}
N1deta2 = 1./4.*(1.0-eta1)*(-1.0);  % N^1_{,eta2}
N2deta1 = 1./4.*(+1.0)*(1.0-eta2);  
N2deta2 = 1./4.*(1.0+eta1)*(-1.0);  
N3deta1 = 1./4.*(+1.0)*(1.0+eta2);
N3deta2 = 1./4.*(1.0+eta1)*(+1.0);
N4deta1 = 1./4.*(-1.0)*(1.0+eta2);
N4deta2 = 1./4.*(1.0-eta1)*(+1.0);

Ndeta = [N1deta1,N2deta1,N3deta1,N4deta1;...
         N1deta2,N2deta2,N3deta2,N4deta2];
     
%% Jacobi Matrix of element e, evaluated at GP (eta)
% J = [dx1/deta1 dx2/deta1;
%      dx1/deta2 dx2/deta2]
xbar_var = [xbar([1:2:8]),xbar([2:2:8])];
J = Ndeta*xbar_var; 
Jinv = J\eye(length(J));
deta1dx1 = Jinv(1,1);
deta1dx2 = Jinv(1,2);
deta2dx2 = Jinv(2,2);
deta2dx1 = Jinv(2,1);

%% compute B-operator

N1x1 = N1deta1*deta1dx1 + N1deta2*deta2dx1; % N^1_{,x1}
N1x2 = N1deta1*deta1dx2 + N1deta2*deta2dx2; % N^1_{,x2}
N2x1 = N2deta1*deta1dx1 + N2deta2*deta2dx1;  
N2x2 = N2deta1*deta1dx2 + N2deta2*deta2dx2;
N3x1 = N3deta1*deta1dx1 + N3deta2*deta2dx1;
N3x2 = N3deta1*deta1dx2 + N3deta2*deta2dx2;
N4x1 = N4deta1*deta1dx1 + N4deta2*deta2dx1;
N4x2 = N4deta1*deta1dx2 + N4deta2*deta2dx2;


B = [N1x1,0,N2x1,0,N3x1,0,N4x1,0;...
     0,N1x2,0,N2x2,0,N3x2,0,N4x2;...
     N1x2,N1x1,N2x2,N2x1,N3x2,N3x1,N4x2,N4x1];
end

