%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Solar/Solar_Functions/');
fprintf('All paths added...\n\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data paths
project_path   = '/Volumes/Solar_3/SecondarySolar_4';
recording_name = 'SS_NOM_S15_H3_A000_B000_PIVMID_UH';
processing     = 'StereoPIV_MPd(2x24x24_50%ov)';
inpt_name      = recording_name;

% Image paths
piv_path = fullfile(project_path, recording_name, processing);

% Save paths
results_path = '/Volumes/Solar_3/Test_Processing/';
mtlb_file    = strcat(results_path, 'data'   , '/', inpt_name, '_DATA.mat');
mean_file    = strcat(results_path, 'means'  , '/', inpt_name, '_MEANS.mat');
figure_file  = strcat(results_path, 'figures', '/', inpt_name);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAVIS TO MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(mtlb_file, 'file')
    fprintf('* Loading DATA from File\n')
    data = load(mtlb_file);
    data = data.output;
else
    data = vector2matlab(piv_path, mtlb_file);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB DATA TO ENSEMBLE/PHASE MEANS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(mean_file, 'file')
     fprintf('* Loading MEANS from File\n')
     means = load(mean_file); 
     means = means.output;
else
     means = data2means(mean_file, data);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
% Coordinates
X = means.X.';
Y = means.Y.';

% Means
U = means.u;
V = means.v;
W = means.w;

% Normal Stresses
uu = means.uu;
vv = means.vv;
ww = means.ww;

% Shear Stresses
uv = means.uv;
uw = means.uw;
vw = means.vw;

%%% MASK
U(U == 0) = nan;
V(U == 0) = nan;
W(U == 0) = nan;


% Means Plots

ax = figure();
t  = tiledlayout(1,3);
sgtitle(inpt_name, 'interpreter', 'none')

nexttile()
colormap jet
contourf(X, Y, U, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('u')

nexttile()
colormap jet
contourf(X, Y, V, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('v')

nexttile()
colormap jet
contourf(X, Y, W, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('w')

%% Stresses Plots

ax = figure();
t  = tiledlayout(2,3);
sgtitle(inpt_name, 'interpreter', 'none')

% Normal Stresses
nexttile()
colormap jet
contourf(X, Y, uu, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('uu')

nexttile()
colormap jet
contourf(X, Y, vv, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('vv')

nexttile()
colormap jet
contourf(X, Y, ww, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('ww')


% Shear Stresses
nexttile()
colormap jet
contourf(X, Y, uv, 100, 'linestyle', 'none')
axis equal
xlim([-100,100])
ylim([-100,100])
colorbar()
title('uv')

nexttile()
colormap jet
contourf(X, Y, uw, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('uw')

nexttile()
colormap jet
contourf(X, Y, vw, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('vw')


%%



figure()
contourf(X, Y, means.w(:,:,1), 500, 'linestyle', 'none')
axis equal
















