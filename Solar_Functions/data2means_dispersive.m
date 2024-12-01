% This function converts the Matlab vector data from PIV to Means 
% and Reynolds Stresses for further processing.
% Ondrej Fercak, Zein Sadek, 3/21/2022
%Abdelhalim Abdeldayem, added dispersive variables 9/4/2024

% out_path:     Folder where new struct file will be saved.
% out_name:     Name of new struct file.
% inst_struct:  Matlab data struct as input for calculations.

function output = data2means_dispersive(out_path, inst_struct)

    % Check if Save Folder Exists. [if not, create]
    if exist(out_path, 'file')
        fprintf('<data2means> *Save Folder was Previously Created. \n')
    else
        fprintf('<data2means> *Creating New Save Floder. \n')
    end
   
    D = inst_struct.D;
    % fprintf('D Check = %d! \n', D)

    % Extract Instantaneous Velocities from Struct.
    inst_u  = inst_struct.U;
    inst_v  = inst_struct.V;
    inst_w  = inst_struct.W;

    % x       = inst_struct.X;
    % y       = inst_struct.Y;
    inst_u(inst_u==0)=NaN;
    inst_v(inst_v==0)=NaN;
    inst_w(inst_w==0)=NaN;

    % Calculate Velocity Means
    fprintf('Computing Ensemble Average\n')
    output.u = mean(inst_u, 3, 'omitnan');
    output.v = mean(inst_v, 3, 'omitnan');
    output.w = mean(inst_w, 3, 'omitnan');

    % Calculate Velocity profiles, average over spanwise direction (z)
    fprintf('Computing Spanwise Average\n')
    output.uz=mean(output.u,2,'omitnan');
    output.vz=mean(output.v,2,'omitnan');
    output.wz=mean(output.w,2,'omitnan');
    
    % Calculate Dispersive Velocity 
    fprintf('Computing Dispersive Velocity\n')
    output.ud=output.u - output.uz;
    output.vd=output.v - output.vz;
    output.wd=output.w - output.wz;
    
    % Create Reynolds Stress Objects
    uu_p = zeros(size(inst_u));
    vv_p = zeros(size(inst_u));
    ww_p = zeros(size(inst_u));
    
    uv_p = zeros(size(inst_u));
    uw_p = zeros(size(inst_u));
    vw_p = zeros(size(inst_u));

    % Loop Through Each Frame in Struct. 
    %{ 
        [Halim]: no need for this loop. Consider doing the following
         u_p = inst_u - output.u; %subtracting a 2d matrix from each frame of the 3d matrix (matlab is smart)
         output.uu=mean(u_p.*u_p,3,'omitnan')
    %}
    % output.u; uu_p
    fprintf('\n<data2means> PROGRESS: ');
    for frame_number = 1:D
        
        % Print Progress.
        progressbarText(frame_number/D);
        
        % Instantaneous Fluctuations.
        u_pi = inst_u(:, :, frame_number) - output.u;
        v_pi = inst_v(:, :, frame_number) - output.v;
        w_pi = inst_w(:, :, frame_number) - output.w;

        % Instantaneous Stresses.
        uu_pi = u_pi.*u_pi;
        vv_pi = v_pi.*v_pi;
        ww_pi = w_pi.*w_pi;

        uv_pi = u_pi.*v_pi;
        uw_pi = u_pi.*w_pi;
        vw_pi = v_pi.*w_pi;

        % Array of Mean Stresses.
        uu_p(:, :, frame_number) = uu_pi;
        vv_p(:, :, frame_number) = vv_pi;
        ww_p(:, :, frame_number) = ww_pi;

        uv_p(:, :, frame_number) = uv_pi;
        uw_p(:, :, frame_number) = uw_pi;
        vw_p(:, :, frame_number) = vw_pi;

    end

    % Turbulent Mean Stresses.
    output.uu = mean(uu_p, 3, 'omitnan');
    output.vv = mean(vv_p, 3, 'omitnan');
    output.ww = mean(ww_p, 3, 'omitnan');

    output.uv = mean(uv_p, 3, 'omitnan');
    output.uw = mean(uw_p, 3, 'omitnan');
    output.vw = mean(vw_p, 3, 'omitnan');
    
    %Dispersive stresses
    
    % Normal Dispersive stresses
    output.uud=output.ud.*output.ud;
    output.vvd=output.vd.*output.vd;
	output.wwd=output.wd.*output.wd;

	output.uvd=output.ud.*output.vd;
    output.uwd=output.ud.*output.wd;
	output.vwd=output.vd.*output.wd;
    
    output.X = inst_struct.X;
    output.Y = inst_struct.Y;
    output.D = D;
    
    % Save Matlab File.
    fprintf('\n<data2means> Saving Data to File... \n');
    save(out_path, 'output');
    clc; fprintf('\n<data2means> Data Save Complete \n')
end
