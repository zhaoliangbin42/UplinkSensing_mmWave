%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Function for detecting the Line-of-Sight (LoS) path in mmWave communications
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y_AoA] = Func_Detect_LoS(Y_AoA, sig_param, FFT_length, is_plot)
    % FUNC_DETECT_LOS - Detects the Line-of-Sight (LoS) path from estimated AoAs
    % Inputs: 
    %   Y_AoA: Struct containing received signal data
    %       - Y: Received signal [Nr x Snap x K] (complex matrix)
    %       - AoA_est: Estimated AoAs [L x 1] (degrees)
    %   sig_param: Structure containing signal parameters
    %       - Snap_AoA: Number of snapshots for AoA estimation (integer)
    %       - Nr: Number of receive antennas (integer)
    %       - K: Number of subcarriers (integer)
    %       - Ls: Number of static paths (integer)
    %       - Ld: Number of dynamic paths (integer)
    %       - W: Combining matrix of the angle codebook (complex matrix)
    %       - AoA: Actual AoAs [L x 1] (degrees)
    %   FFT_length: Length of the FFT for spatial processing (integer)
    %   is_plot: Flag to enable plotting and verbose output (boolean)
    % Outputs: 
    %   Y_AoA: Updated struct with added field:
    %       - LoS_AoA_est: Estimated LoS AoA (degrees)
    % Example: 
    %   [Y_AoA] = Func_Detect_LoS(Y_AoA, sig_param, 1024, true)

    % Extract parameters from input structures
    Snap        = sig_param.Snap_AoA;
    Nr          = sig_param.Nr;
    K           = sig_param.K;
    Ls          = sig_param.Ls;
    Ld          = sig_param.Ld;
    W           = sig_param.W;
    AoA         = sig_param.AoA;
    b           = sig_param.b;

    Y           = Y_AoA.Y;
    AoA_est     = Y_AoA.AoA_est;
    L           = Ls + Ld;  % Total number of paths

    % Extract first snapshot and subcarrier to recover channel information
    first_snap  = squeeze(Y(1, :, 1));

    % Recover the channel state information vector by applying IDFT
    % NOTE: W is the IDFT matrix that transforms from angle domain to spatial domain
    h           = W * first_snap(:);

    % Zero padding to increase angle resolution for FFT processing
    h_padded    = [h; zeros(FFT_length - length(h), 1)]; 

    % Generate array response pattern using FFT
    H           = Func_Array_Response(h_padded, FFT_length, is_plot);
    if is_plot
        hold on;
        % Plot static path (blue) and dynamic paths (red)
        stem(AoA(1), abs(b(1)), 'bo', 'MarkerSize', 10, 'LineWidth', 1.5);
        stem(AoA(2:end), abs(b(2:end)), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
        hold off;
    end

    % Find the peak in the array response pattern (strongest path direction)
    [~, idx]    = max(abs(H));
    
    % Create angle codebook mapping FFT bins to physical angles
    % NOTE: u = sin(¦È)/2 in the spatial frequency domain
    u_codebook      = (-FFT_length/2 : FFT_length/2 - 1) / FFT_length;
    angle_codebook  = asin(u_codebook * 2) * 180 / pi;
    
    % Get the beam angle corresponding to the peak response
    LoS_beam_angle  = angle_codebook(idx);

    % Find the estimated AoA closest to the beam angle (LoS path)
    [~, idx]        = min(abs(AoA_est - LoS_beam_angle));
    LoS_AoA_est     = AoA_est(idx);

    % Update output structure with estimated LoS AoA
    Y_AoA.LoS_AoA_est = LoS_AoA_est;
    
    % Display results if plotting is enabled
    if is_plot
        disp("---------------------- LoS Detection ----------------------");
        disp(sprintf("The estimated LoS AoA is %f", LoS_AoA_est));
        
        % Check if LoS detection was successful (within 1 degree tolerance)
        isLoSDetected = Func_Compare_Sets(LoS_AoA_est, AoA(1), 1);
        
        if isLoSDetected
            disp('Conclusion: LoS path detected successfully.');
        else
            disp('Conclusion: LoS path detection failed.');
        end
        % disp("-----------------------------------------------------------");
    end
end