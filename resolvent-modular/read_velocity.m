%% Sheel Nidhan
%  Date - 29th January 2020

clear; clc;
%% Read grid

fid = fopen('/home/sheel/Work/projects/spod_re5e4/grid/frinf/x1_grid.in');  %% Reading the radial grid

D = cell2mat(textscan(fid, '%f%f', 'headerlines', 1));

r = D(1:end-9,2);

for i = 1:size(r,1)-2
    rc(i,1) = 0.5*(r(i+1,1) + r(i,1));  % Centered the grid faces to grid centers
end


%% Read velocity 

load('./ustreamwise/mean_velocity_x_D_100.mat');
load('./ustreamwise/mean_velocity_x_D_30.mat');
