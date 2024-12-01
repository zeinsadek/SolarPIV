%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx/');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Solar/Solar_Functions/');
fprintf('All Paths Imported...\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data paths
clc;
project_path   = '/Volumes/Solar_APS/CleanDriveTest';
recording_name = 'SS_NOM_S15_H3_A045_B045_PIVSPN_UH';
processing     = 'StereoPIV_MPd(2x24x24_50%ov)';
inpt_name      = strcat('D_', recording_name, '_REMAP');

% Image paths
piv_path = fullfile(project_path, recording_name, processing);

% Save paths
results_path = '/Volumes/Solar_APS/Results/';
mtlb_file    = strcat(results_path, 'data'   , '/', inpt_name, '_DATA.mat');
mean_file    = strcat(results_path, 'means'  , '/', inpt_name, '_MEANS.mat');
figure_file  = strcat(results_path, 'figures', '/', inpt_name);

% Make specific folder for figures of an experiment
% if ~exist(figure_file, 'dir')
%     mkdir(figure_file)
% end

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
     means = data2means_dispersive(mean_file, data);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
% Coordinates
X = means.X;
Y = means.Y;

% Means
U = means.u;
V = means.v;
W = means.w;

% Velocity profiles, average over spanwise direction (y)
Uz=means.uz;
Vz=means.vz;
Wz=means.wz;

% Dispersive velcoity
Ud=means.ud;
Vd=means.vd;
Wd=means.wd;

Ud(Ud==0)=NaN;
Vd(Vd==0)=NaN;
Wd(Wd==0)=NaN;

% Normal Turbuluent Stresses
uu = means.uu;
vv = means.vv;
ww = means.ww;

uu(U==0)=NaN;
vv(U==0)=NaN;
ww(U==0)=NaN;

% Shear turbulent Stresses
uv = means.uv;
uw = means.uw;
vw = means.vw;

uv(U==0)=NaN;
uw(U==0)=NaN;
vw(U==0)=NaN;

% Normal Dispersive stresses
uud=means.uud;
vvd=means.vvd;
wwd=means.wwd;

uud(U==0)=NaN;
vvd(U==0)=NaN;
wwd(U==0)=NaN;

% Shear Dispersive stresses
uvd=means.uvd;
uwd=means.uwd;
vwd=means.vwd;

uvd(U==0)=NaN;
uwd(U==0)=NaN;
vwd(U==0)=NaN;

V(U==0)=NaN;
W(U==0)=NaN;
U(U==0)=NaN;
%% Means Plots

kernel = 20;

ax = figure();
t  = tiledlayout(1,3);
sgtitle(inpt_name, 'interpreter', 'none')

nexttile()
colormap jet
hold on
contourf(X, Y, juliaan_smooth(U, kernel), 100, 'linestyle', 'none')
hold off
axis equal
clim([0, 4])
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('u')

nexttile()
colormap jet
contourf(X, Y, juliaan_smooth(V, kernel), 100, 'linestyle', 'none')
axis equal
clim([-0.5, 0.5])
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('v')

nexttile()
colormap jet
contourf(X, Y, juliaan_smooth(W, kernel), 100, 'linestyle', 'none')
axis equal
clim([-2.5 0.5])
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('w')




%% Dispersive Velocity Plots

kernel = 25;
Wd = juliaan_smooth(Wd, kernel);
Vd = juliaan_smooth(Vd, kernel);

vn=5;
ax = figure();
%t  = tiledlayout(1,3);
sgtitle(inpt_name, 'interpreter', 'none')

nexttile()
colormap jet
contourf(X, Y, Ud, 500, 'linestyle', 'none')
colormap('coolwarm')
axis equal
hold on
% xlim([-100,100])
% ylim([-100,100])
colorbar()
title('ud')
clim([-1, 1]);
quiver(X(1:vn:end,1:vn:end),Y(1:vn:end,1:vn:end),Wd(1:vn:end,1:vn:end),Vd(1:vn:end,1:vn:end), 2.5,'color',[0 0 0]);

% nexttile()
% colormap jet
% contourf(X, Y, uvd, 500, 'linestyle', 'none')
% colormap('coolwarm')
% axis equal
% hold on
% % xlim([-100,100])
% % ylim([-100,100])
% colorbar()
% title('ud')
% clim([-1E-3, 1E-3]);
% quiver(X(1:vn:end,1:vn:end),Y(1:vn:end,1:vn:end),Vd(1:vn:end,1:vn:end),Wd(1:vn:end,1:vn:end), 2.5,'color',[0 0 0]);

%nexttile()
% colormap jet
% contourf(X, Y, Vd, 100, 'linestyle', 'none')
% axis equal
% hold on
% xlim([-100,100])
% ylim([-100,100])
% colorbar()
% title('vd')
%     quiver(X(1:vn:end,1:vn:end),Y(1:vn:end,1:vn:end),Wd(1:vn:end,1:vn:end),Ud(1:vn:end,1:vn:end),10,'color',[0 0 0]);
% 
% nexttile()
% colormap jet
% contourf(X, Y, Wd, 100, 'linestyle', 'none')
% axis equal
% hold on
% xlim([-100,100])
% ylim([-100,100])
% colorbar()
% title('wd')
%     quiver(X(1:vn:end,1:vn:end),Y(1:vn:end,1:vn:end),Ud(1:vn:end,1:vn:end),Vd(1:vn:end,1:vn:end),10,'color',[0 0 0]);

%% Velocity profiles

ax = figure();
%t  = tiledlayout(1,3);
sgtitle(inpt_name, 'interpreter', 'none')

UUU=mean(U,1,'omitnan');

nexttile()
hold on
plot(Uz,Y(1,:))
plot(UUU,Y(1,:))
title('Dispersive u')


for i=1:2:200
    plot(U(i,:),Y(1,:))
end

% nexttile()
% plot(Vz,Y)
% title('Dispersive v')
% 
% nexttile()
% plot(Wz,Y)
% title('Dispersive w')

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
%[Halim] Consider using redblue colormap for shear stresses. You need to place the redblue.m file in the folder
nexttile()
colormap jet
contourf(X, Y, uv, 100, 'linestyle', 'none')
axis equal
% xlim([-100,100])
% ylim([-100,100])
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


%% Dispersive Stresses Plots

ax = figure();
t  = tiledlayout(2,3);
sgtitle(inpt_name, 'interpreter', 'none')

% Dispersive Normal Stresses
nexttile()
colormap jet 
contourf(X, Y, uud, 100, 'linestyle', 'none')
axis equal
xlim([-100,100])
ylim([-100,100])
colorbar()
title('Dispersive uu')
colormap('coolwarm')

nexttile()
colormap jet
contourf(X, Y, vvd, 100, 'linestyle', 'none')
axis equal
xlim([-100,100])
ylim([-100,100])
colorbar()
title('Dispersive vv')
colormap('coolwarm')

nexttile()
colormap jet
contourf(X, Y, wwd, 100, 'linestyle', 'none')
axis equal
xlim([-100,100])
ylim([-100,100])
colorbar()
title('Dispersive ww')
colormap('coolwarm')


% Dispersive Shear Stresses
%[Halim] Consider using redblue colormap for shear stresses. You need to place the redblue.m file in the folder
nexttile()
colormap jet 
contourf(X, Y, uvd, 100, 'linestyle', 'none')
axis equal
xlim([-100,100])
ylim([-100,100])
colorbar()
title('Dispersive uv')
colormap('coolwarm')

nexttile()
colormap jet
contourf(X, Y, uwd, 100, 'linestyle', 'none')
axis equal
xlim([-100,100])
ylim([-100,100])
colorbar()
title('Dispersive uw')
colormap('coolwarm')

nexttile()
colormap jet
contourf(X, Y, vwd, 100, 'linestyle', 'none')
axis equal
xlim([-100,100])
ylim([-100,100])
colorbar()
title('Dispersive vw')
colormap('coolwarm')



















