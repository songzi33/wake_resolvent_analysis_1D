%% Name - Sheel Nidhan
%  Transformation based on Lesshafft and Huerre 2007 eqn. A13a
%  Just in testing phase, consolidate later

close all;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');  set(groot, 'defaultTextInterpreter','latex'); 
set(groot, 'defaultFigureRenderer','painters')
set(groot, 'defaultFigureColor',[1 1 1])

%% Checking accuracy of Chebysev differentiation (Matching)


f = x(half).^2;

figure;
plot(x(half), f, 'r-', 'Linewidth', 2);
hold on;
plot(x(half), D1(half,half)*f(half), 'k-', 'Linewidth', 2);  %% Looks accurate
% plot(x(half), D2(half, half)*f, 'b-', 'Linewidth', 2);  %% Looks accurate

%% Transform the x to 0,rmax as given in Lesshafft and Huerre 2007 (Matching)

rc = 3; rmax = 10;

r_wake = rc*(1-x)./(1-x.^2 + 2*rc/rmax);

f = r_wake.^2;

figure;
plot(r_wake, f, 'r-', 'Linewidth', 2);
hold on;

% transformation vector needed

drdeta = (rc*(2*x - 1- x.^2 - 2*rc/rmax)./((1-x.^2 + 2*rc/rmax).^2));
detadr = 1./drdeta;
detadr_diag = diag(detadr);
D1_transformed = detadr_diag*D1;
plot(r_wake, D1_transformed*f, 'k-', 'Linewidth', 2);
plot(r_wake, 2*r_wake, 'bs');

%% Testing the second derivative based on Lesshaft and Huerre 2007 (Matching)

detadrdeta = ((1-x.^2 + 2*rc/rmax).^2*rc.*(2-2.*x) - rc*(2*x-1-x.^2-2*rc/rmax).*2.*(1-x.^2+2*rc/rmax).*(-2*x))./((1-x.^2 + 2*rc/rmax).^4);
detadrdeta_diag = diag(detadrdeta);
% d2etadr2 = 2*((1-x.^2 + 2*rc/rmax).^3).*((1-x).*(1-x.^2+2*rc/rmax) + 2*x.*(2*x-1-x.^2-2*rc/rmax))./(rc.^2.*(2*x-1-x.^2-2*rc/rmax).^3);
% d2etadr2_diag = diag(d2etadr2);

D2_transformed = -(detadrdeta_diag)*(detadr_diag.^3)*D1 + (detadr_diag.^2)*D2;

hold on;
plot(r_wake, D2_transformed*f, 'k-', 'Linewidth', 2);
plot(r_wake, 2, 'ro');

ylim([0 2.1]);
    