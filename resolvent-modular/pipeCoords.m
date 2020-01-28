function [r, dr, D1E, D1O, D2E, D2O] = pipeCoords( n, N, rc2, rmax)
% Function to set up the coordinate system, differentiation and integration
% matrices for the pipe geometry.

% Written by Mitul Luhar, 02/06/2013

% r: radial coordinate (0,1]
% dr: integration weights
% D1E,D1O:  Even and odd first difference matrices
% D2E,D2O:  Even and odd second difference matrices

[x,DM] = chebdif(N,2); % x:[-1 1]
half   = 1:N;            % indices corresponding to x:(0,1]

r = 0.001 + rc2*(1-x)./(1-x.^2 + 2*rc2/rmax); % r corresponding to the wake crossection
r = r(1:N);
% Transformation matrices needed for D1 and D2
drdeta = (rc2*(2*x - 1- x.^2 - 2*rc2/rmax)./((1-x.^2 + 2*rc2/rmax).^2));
detadr = 1./drdeta(1:N);
detadr_diag = diag(detadr);

detadrdeta = ((1-x.^2 + 2*rc2/rmax).^2*rc2.*(2-2.*x) - rc2*(2*x-1-x.^2-2*rc2/rmax).*2.*(1-x.^2+2*rc2/rmax).*(-2*x))./((1-x.^2 + 2*rc2/rmax).^4);
detadrdeta_diag = diag(detadrdeta(1:N));

D1 = DM(1:N,1:N,1);          % First differential
D2 = DM(1:N,1:N,2);          % Second differential

%Split D1 and D2 into odd and even components (cf. Meseguer Trefethen 2003)
%For even (odd) n, u is even (odd) over r:[-1 1] 
%For even (odd) n, v and w are odd (even) over r:[-1 1]
s = (-1)^mod(n,2);
D1E = detadr_diag*(D1 + s*D1);
D1O = detadr_diag*(D1 - s*D1);
D2E = ((detadr_diag.^2)*D2 -(detadrdeta_diag)*(detadr_diag.^3)*D1) +  ...
       s*(((detadr_diag).^2*D2 -(detadrdeta_diag)*(detadr_diag.^3)*D1));
D2O = ((detadr_diag.^2)*D2 -(detadrdeta_diag)*(detadr_diag.^3)*D1) -  ...
       s*(((detadr_diag.^2)*D2 -(detadrdeta_diag)*(detadr_diag.^3)*D1));
%D1E = D1(half,half) + s*D1(half,2*N+1-half);
%D1O = D1(half,half) - s*D1(half,2*N+1-half);
%D2E = D2(half,half) + s*D2(half,2*N+1-half);
%D2O = D2(half,half) - s*D2(half,2*N+1-half);

% Find integration weights %change later
[~,dr] = clencurt(2*N-1);
dr   = (dr(half))';

end

