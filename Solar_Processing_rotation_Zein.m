% Solar Perspective Projected PIV

clc
clear all
close all

addpath("/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Solar/Solar_Functions")
addpath("/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Solar/Solar_Functions/AxelRot")

%% Setup

% Case name: B45 is beta 45 deg
means = load("/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Solar/data/D_SS_NOM_S15_H3_A045_B045_PIVSPN_UH_REMAP_MEANS.mat");
means = means.output;

% Parameters
U_inf = 4;
U_inf_sq = U_inf^2;

% Means
U = means.u./U_inf;
V = means.w./U_inf;
W = means.v./U_inf;

% Dispersive velcoity
Ud = means.ud./U_inf;
Vd = means.vd./U_inf;
Wd = means.wd./U_inf;

%% Coordinates (2D mesh)

% PIV coordinates (2D)
ZZ = means.X;
YY = means.Y;

% PIV coordinates (1D)
z = ZZ(1,:);
y = YY(:,1);

% Grid spacing in y and z
dyz = abs(y(1)-y(2)); 

% Building the x coodrinate 
x = 0:dyz:170*dyz;

% Building the 3D mesh
[ZZZ,YYY,XXX] = meshgrid(z,y,x); 

% Rotating to the new coordinates
Udold = reshape(Ud, 1, []);
Vdold = reshape(Vd, 1, []);
Wdold = reshape(Wd, 1, []);

XXXold = reshape(XXX, 1, []);
YYYold = reshape(YYY, 1, []);
ZZZold = reshape(ZZZ, 1, []);

%% Coordinates (3D mesh)

[temp,~,~] = AxelRot([Udold;Vdold;Wdold], -45, [0 1 0]', []);

Udnew = reshape(temp(1,:), [length(y), length(z)]);
Vdnew = reshape(temp(2,:), [length(y), length(z)]);
Wdnew = reshape(temp(3,:), [length(y), length(z)]);

[temp,~,~] = AxelRot([XXXold;YYYold;ZZZold], -45, [0 1 0]', []);

XXXnew = reshape(temp(1,:), [length(y), length(z), length(x)]);
YYYnew = reshape(temp(2,:), [length(y), length(z), length(x)]);
ZZZnew = reshape(temp(3,:), [length(y), length(z), length(x)]);

%% Rebuilding velocity vectors

% Building the 3d velocity matrises
Ud3d = 0*XXX;
Vd3d = 0*XXX;
Wd3d = 0*XXX;

% Building the 3D velocity matrices
for i = 1:length(x)

    Ud3d(:,1:end-(i-1),i) = Udnew(:,i:end);
    Ud3d(:,end-(i-1)+1:end,i) = Udnew(:,1:i-1);

    Vd3d(:,1:end-(i-1),i) = Vdnew(:,i:end);
    Vd3d(:,end-(i-1)+1:end,i) = Vdnew(:,1:i-1);

    Wd3d(:,1:end-(i-1),i) = Wdnew(:,i:end);
    Wd3d(:,end-(i-1)+1:end,i) = Wdnew(:,1:i-1);
end

% Extracting the velocities of the rotated plane 
Udrot = zeros(length(y), length(x));
Vdrot = zeros(length(y), length(x));
Wdrot = zeros(length(y), length(x));

% Extracting the tilted plane from the rotated matrix 
for i=1:length(x)

    Udrot(:,i) = Ud3d(:,i,i);
    Vdrot(:,i) = Vd3d(:,i,i);
    Wdrot(:,i) = Wd3d(:,i,i);

    % XXXrot(:,i)=XXXnew(:,i,i);
    % YYYrot(:,i)=YYYnew(:,i,i);
    % ZZZrot(:,i)=ZZZnew(:,i,i);
     YYrot(:,i) = YYYnew(:,i,i);
     ZZrot(:,i) = ZZZnew(:,i,i);
end


%% Dispersive Velocity Plots Rotated

% kernel = 25;
% Wd = juliaan_smooth(Wd, kernel);
% Vd = juliaan_smooth(Vd, kernel);

Reference_Velocity = 0.05;
Scale = double(max(max(sqrt(Vd.^2 + Wd.^2 )))/Reference_Velocity);

figure();
vn = 2;
l_Du = 0.05;

nexttile()
contourf(ZZrot, YYrot, Udrot, 500, 'linestyle', 'none')
% colormap('redblue')
axis equal
hold on
xlim([-150,150])
ylim([-100,100])
% xlim([-x_crop, x_crop])
% ylim([0, y_crop])
colorbar()
title('$<\overline{u}> \,/ U_{\infty}^2$', 'Interpreter', 'latex','FontSize', 20)
clim([-l_Du l_Du]);
quiver(ZZrot(1:vn:end,1:vn:end),YYrot(1:vn:end,1:vn:end),Wdrot(1:vn:end,1:vn:end),Vdrot(1:vn:end,1:vn:end),Scale,'color',[0 0 0]);
ylabel('$y/S_y$', 'Interpreter', 'latex')
xlabel('$z/S_y$', 'Interpreter', 'latex')
xlabel('z','FontSize', 20)
ylabel('y','FontSize', 20)

% print(gcf,'Dew.svg','-dsvg','-r1200');  
%exportgraphics(gcf,strcat(cs,'_rot.png'),'Resolution',3000) 

%out_path = '/Users/katietaylor/Library/Mobile Documents/com~apple~CloudDocs/Research/PIV/Solar/Results/figures/Velocity';
%exportgraphics(ax, fullfile(out_path, 'Dispersive Velocity_B0.png'), 'Resolution', 200);

%{
%% plotting in 3D

figure ()

for i=1:10

    p=patch(isosurface(ZZZ,YYY,XXX,Ud3d,-0.5*i/1000));
    f = i/10; % fraction number
    p.FaceColor = [1 f 1-f];
    p.EdgeColor = 'none';
    p.FaceAlpha = 1-f;
    axis tight
    camlight 
    lighting gouraud
end

daspect([1 1 0.1])

view(0,91)
camroll(180)

for i=1:36
   camorbit(1,1,'camera')
   drawnow
   pause(0.25)
end
% view(100,15)
% camroll(90)
%}

%% Dispersive Velocity Plots

% kernel = 25;
% Wd = juliaan_smooth(Wd, kernel);
% Vd = juliaan_smooth(Vd, kernel);

Reference_Velocity=0.05;
Scale=double(max(max(sqrt(Vd.^2 + Wd.^2 )))/Reference_Velocity);

ax = figure();
%t  = tiledlayout(1,3);
% sgtitle(inpt_name, 'interpreter', 'none')
%l_Du = 0.2;
l_Du = 0.05;
vn=2;

nexttile()
contourf(ZZ, YY, -Ud, 500, 'linestyle', 'none')
colormap('redblue')
axis equal
hold on
xlim([-200,150])
ylim([-100,100])
% xlim([-x_crop, x_crop])
% ylim([0, y_crop])
colorbar()
title('$<\overline{u}> \,/ U_{\infty}^2$', 'Interpreter', 'latex','FontSize', 20)
clim([-l_Du l_Du]);
quiver(ZZ(1:vn:end,1:vn:end),YY(1:vn:end,1:vn:end),Wd(1:vn:end,1:vn:end),Vd(1:vn:end,1:vn:end),Scale,'color',[0 0 0]);
ylabel('$y/S_y$', 'Interpreter', 'latex')
xlabel('$z/S_y$', 'Interpreter', 'latex')
xlabel('z','FontSize', 20)
ylabel('y','FontSize', 20)
%patch([-1.2 -1.2 1.02 1.02],[-0.02 0 0 -0.02 ],'k')
%patch([-200 -200 150 150],[-105 -100 -100 -105 ],'k')

%exportgraphics(gcf,strcat(cs,'.png'),'Resolution',3000) 
print(gcf,'Nor.svg','-dsvg','-r1200');  

%out_path = '/Users/katietaylor/Library/Mobile Documents/com~apple~CloudDocs/Research/PIV/Solar/Results/figures/Velocity';
%exportgraphics(ax, fullfile(out_path, 'Dispersive Velocity_B0.png'), 'Resolution', 200);
