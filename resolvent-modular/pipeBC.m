function [H] = pipeBC(LHS,RHS,N,r,m)
% Function to impose boundary conditions for Navier-Stokes Resolvent
% Written by Mitul Luhar on 02/06/2013

% After Fourier decomposition, the Navier-Stokes eqns are:
% (-1i*om*M - L)x = Mf , or LHS x = RHS f
% where x = [u;v;w;p] and f = [fu;fv;fw;-]

% varargin specifies the boundary condition
% varargin = {} - no slip
% varargin = {yPD,AD}: opposition control, detection at yPD (plus units)
% and amplitude AD, such that v(yPD) = -AD*v(0)

% yP: coordinates (in plus units)

% Implement boundary conditions by changing rows 1, N+1, 2*N+1 in LHS/RHS. 
% These correspond to the x,r,theta momentum balance at r = 1 (y = 0). 

%% Boundary condition at r = 0
if m == 1 
    % Boundary condition for m = 1

    index_m1 = [1 3*N+1];
    LHS(index_m1,:) = 0;  %ux=p=0 at r=0
    RHS(index_m1,:) = 0;
        LHS(1, 1) = 1; LHS(3*N+1, 3*N+1) = 1;

    LHS(2*N+1,:) = 0;
    LHS(2*N+1,2*N+1) = sqrt(-1); LHS(2*N+1,N+1) = 1;   % utheta + iur = 0
    RHS(2*N+1,:) = 0;

    LHS(N+1, :) = 0; % ur' = 0
    RHS(N+1, :) = 0;  

    LHS(N+1,N+1+1) = -1/(r(3)-r(2));
    LHS(N+1,N+1+2) =  1/(r(3)-r(2));

elseif m > 1

    % Boundary condition at axis for |m| > 1
    % Set rows [1,N+1,2*N+1, 3*N+1] in LHS = 0.  To be modified below to reflect BC
    LHS(1:N:4*N,:) = 0;
    % Set rows [1,N+1,2*N+1, 3*N+1] in RHS = 0.
    RHS(1:N:4*N,:) = 0;

    LHS(1,1) = 1;
    LHS(N+1,N+1)=1;
    LHS(2*N+1,2*N+1) = 1;
    LHS(3*N+1,3*N+1) = 1;

end

%% Boundary condition at r = 10D

% Set rows [1,N+1,2*N+1] in LHS = 0.  To be modified below to reflect BC
index_set = [N 2*N 3*N 4*N];
LHS(index_set,:) = 0;
% Set rows [1,N+1,2*N+1] in RHS = 0.
RHS(index_set,:) = 0;

LHS(N,N)     = 1;
LHS(2*N,2*N) = 1;
LHS(3*N,3*N) = 1;
LHS(4*N,4*N) = 1;

% if(size(varargin,2)>1)
%     % Opposition control
%     yPD = varargin{1};
%     AD  = varargin{2};
%     ci = find(yP>yPD,1);
%     LHS(N+1,N+1) = 1;
%     LHS(N+1,N+ci)= AD*1;
%     % Row N+1 of mom. balance now represents: v(y=0) + AD*v(y=yD) = 0
% else
%     % No slip
%     LHS(N+1,N+1) = 1;
%     % Row N+1 of mom. balance now represents: v(y=0) = 0
% end

% Compute resolvent which reflects these BCs
% H = pinv(LHS)*RHS;

H = LHS\RHS;

end