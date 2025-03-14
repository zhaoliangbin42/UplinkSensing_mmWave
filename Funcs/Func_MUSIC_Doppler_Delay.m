%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: AoA-based 2D-FFT-MUSIC (AB2FM) using known AoAs
%                 with frequency smoothing
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delay_est, doppler_est] = Func_MUSIC_Doppler_Delay(Y_s, sig_param, grid_num, disp_pixel, is_plot, s_title)
    % FUNC_MUSIC_DOPPLER_DELAY - Performs AB2FM algorithm for Doppler and delay estimation
    % Inputs:  
    %   Y_s: Received signal matrix (Nr*K x Snap_DD)
    %   sig_param: Structure containing signal parameters
    %       - Snap_DD: Number of snapshots
    %       - Nr: Number of receive antennas
    %       - K: Number of subcarriers
    %       - Ls: Number of static paths
    %       - Ld: Number of dynamic paths
    %       - W: Combining matrix
    %       - AoA: Angles of arrival for each path [rad]
    %       - Tm: Symbol duration [s]
    %       - Ts: Sampling time [s]
    %       - fc: Carrier frequency [Hz]
    %       - D_f: Subcarrier spacing [Hz]
    %       - c: Speed of light [m/s]
    %       - Delay_list: List of true delays [s]
    %       - Doppler_list: List of true Doppler shifts [Hz]
    %       - max_distance: Maximum expected distance [m]
    %       - max_velocity: Maximum expected velocity [m/s]
    %   is_plot: Boolean flag for plotting results
    %   s_title: String for plot title
    % Outputs: 
    %   delay_est: Estimated delays [s]
    %   doppler_est: Estimated Doppler shifts [Hz]
    
    % Extract parameters from the structure
    Snap_DD    = sig_param.Snap_DD;
    Nr         = sig_param.Nr;
    K          = sig_param.K;
    Ls         = sig_param.Ls;
    Ld         = sig_param.Ld;
    W          = sig_param.W;
    AoA        = sig_param.AoA;
    Tm         = sig_param.Tm;
    Ts         = sig_param.Ts;
    fc         = sig_param.fc;
    D_f        = sig_param.D_f;
    c          = sig_param.c;
    LoS_Delay  = sig_param.Delay_list(1);  % First delay is LoS reference
    max_distance = sig_param.max_distance;
    max_velocity = sig_param.max_velocity;

    % Total number of paths
    L = Ls + Ld;

    % Frequency smoothing for improved covariance matrix estimation
    % Use half of the subcarriers for smoothing
    sub_K = K / 2;
    
    % Reshape received signal for processing
    H = reshape(reshape(Y_s, Nr, K * Snap_DD), Nr, K, Snap_DD);
    H = reshape(permute(H, [2, 1, 3]), K*Nr, Snap_DD);
    
    % Initialize index array for first smoothing window
    % NOTE: it can be pre-allocated for efficiency, but Nr is typically small
    index = [];
    for i = 0:Nr-1
        index = [index, (i*K+1):(i*K+sub_K)];
    end
    
    % Compute frequency smoothed covariance matrix
    % SEE: Equ. (39)
    R = zeros(Nr*sub_K, Nr*sub_K);
    for k = 1:K-sub_K
        % NOTE: Each smoothing window provides one "snapshot" for covariance estimation
        R = R + H(index, :) * H(index, :)' / Snap_DD;
        index = index + 1;
    end

    % Transfer data from GPU if necessary and compute signal subspace
    R = gather(R);
    [Us_DD, ~] = eigs(R, L, 'largestabs');  % Signal subspace eigenvectors
    
    % Define search grid
    X = grid_num;  % Number of delay points
    Y = grid_num;  % Number of Doppler points
    % SEE: Equ. (52)
    x_axis = linspace(0, max_distance / c, X);            % Delay axis [s]
    y_axis = linspace(-max_velocity, max_velocity, Y) / c * fc / 2;  % Doppler axis [Hz]
    
    UU = Us_DD;  % Signal subspace for MUSIC algorithm

    % Pre-allocate MUSIC spectrum array
    p_MUSIC = zeros(L, X, Y);

    % % Pre-compute steering vectors for Doppler and delay
    % % NOTE: b_f is the candidate vector for Doppler estimation
    % % Not used in this function
    % b_f = (exp(1j * 2 * pi * y_axis.' * Ts .* (0:Nr-1) * Tm).');
    % % NOTE: b_tau is the candidate vector for delays estimation
    % % Not used in this function
    % b_tau = (exp(-1j * 2 * pi * x_axis.' .* (0 + (1:sub_K) * D_f))).';

    % Print header for results
    if is_plot
        s = sprintf('--------------- Estimated Doppler and delay (BF: %s): ---------------', s_title);
        disp(s);
        fprintf('Delay\t\t\t Doppler\n')
    end
    
    % Process each path, implementation of AB2FM algorithm
    for l = 1:L
        p_temp = zeros(X, Y);

        % Get steering vector for the current AoA
        A_dyn = steervec(0.5 * [0:Nr - 1], AoA(l).');
        
        % Apply combining matrix to steering vector
        % SEE: Equ. (51) the left part
        b_w = reshape(repmat(W' * A_dyn, 1, sub_K).', [], 1);

        % Compute MUSIC spectrum using signal subspace
        for i = 1:size(UU, 2)
            % SEE: Equ. (47)
            u = UU(:, i);

            % SEE: Equ. (51)
            u = (u) .* conj(b_w);

            % Reshape into matrix form
            U = reshape(u, sub_K, Nr);

            % NOTE: Use 2D-FFT for efficient computation
            a = fliplr(fftshift(fft2(conj(U), X, Y), 2));
            
            % SEE: The denominator of Equ. (44)
            p_temp = p_temp + abs(a).^2;
        end
        
        % Compute MUSIC pseudospectrum (signal subspace approach)
        % SEE: Equ. (44)
        p_MUSIC(l, :, :) = 1 ./ ((K/2 * Nr) - p_temp);
    end

    % Identify peaks in the 2D MUSIC spectrum
    delay_est_index = zeros(L, 1);
    doppler_est_index = zeros(L, 1);
    
    for i = 1:L
        temp = squeeze(p_MUSIC(i, :, :));
        absTemp = abs(temp);

        % Find global maximum
        [~, maxIndex] = max(absTemp(:));

        % Convert linear index to 2D subscripts
        [delay_est_index(i), doppler_est_index(i)] = ind2sub(size(temp), maxIndex);
        
        if is_plot
            fprintf('%e s\t\t %e Hz\n', x_axis(delay_est_index(i))+LoS_Delay, y_axis(doppler_est_index(i)));
        end
    end
    
    % Convert indices to physical values
    delay_est = x_axis(delay_est_index) + LoS_Delay;
    doppler_est = y_axis(doppler_est_index);

    %% Plot the 2D MUSIC spectrum if requested
    if is_plot
        % Configure plot layout
        row_num = ceil(L/2);
        height_per_row = 0.2;
        figure('units', 'normalized', 'outerposition', [0 0 0.5 row_num*height_per_row]);
        
        % Plot 2D MUSIC spectrum for each path
        for l = 1:L
            subplot(row_num, 2, l);
            
            % Use imresize for efficient and reliable downsampling
            temp_spectrum = squeeze(p_MUSIC(l, :, :));
            p_MUSIC_downsampled = imresize(abs(temp_spectrum), [disp_pixel, disp_pixel], 'nearest');
            
            % Create downsampled coordinate axes for plotting
            x_axis_down = linspace(min(x_axis), max(x_axis), disp_pixel);
            y_axis_down = linspace(min(y_axis), max(y_axis), disp_pixel);
            
            % Plot the spectrum
            mesh(y_axis_down, x_axis_down+LoS_Delay, p_MUSIC_downsampled);
            
            % Mark the peak with a red dot
            hold on
            max_z_value = max(p_MUSIC_downsampled(:));
            plot3(y_axis(doppler_est_index(l)), x_axis(delay_est_index(l))+LoS_Delay, ...
                  max_z_value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            
            % Add text label with peak coordinates
            text(y_axis(doppler_est_index(l)), x_axis(delay_est_index(l))+LoS_Delay, max_z_value, ...
                 sprintf('Delay: %e s\nDoppler: %e Hz', x_axis(delay_est_index(l))+LoS_Delay, ...
                 y_axis(doppler_est_index(l))), 'FontSize', 16, 'FontName', 'Times New Roman');
            
            % Add axis labels and title
            if l < Ls + 1
                s = sprintf('2D MUSIC Spectrum of the static path - %d', l);
            else
                s = sprintf('2D MUSIC Spectrum of the dynamic path - %d', l-Ls);
            end
            Customed_Figure('Doppler(Hz)', 'Delay(s)', s);
            hold off
        end
        
        % Add overall title
        p_title = sprintf('Doppler-delay spectrum, BF: %s', s_title);
        sgtitle(p_title, 'FontSize', 30, 'FontName', 'Times New Roman');
    end
end