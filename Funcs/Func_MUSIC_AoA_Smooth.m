%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Benko
%    Email: lb_zhao_bit_ee@163.com
%
%    Description: Function for AoA estimation using MUSIC algorithm with frequency smoothing
%
%    Tool versions: Matlab 2025a
%    Affiliation: Beijing Institute of Technology
%    Last update: 2025-03-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y_AoA] = Func_MUSIC_AoA_Smooth(Y_AoA, sig_param, grid_num, is_plot)
    % FUNC_MUSIC_AOA_SMOOTH - Angle of Arrival estimation using MUSIC algorithm with frequency smoothing
    % Inputs: 
    %   Y_AoA: Struct containing received signal data
    %       - Y: Received signal [Nr x Snap x K] (complex matrix)
    %   sig_param: Structure containing signal parameters
    %       - Snap_AoA: Number of snapshots for AoA estimation (integer)
    %       - Nr: Number of receive antennas (integer)
    %       - K: Number of subcarriers (integer)
    %       - Ls: Number of static paths (integer)
    %       - Ld: Number of dynamic paths (integer)
    %       - W: Combining matrix of the angle codebook (complex matrix)
    %       - AoA: Actual AoAs [L x 1] (degrees)
    %       - b: Path gains [L x 1] (complex)
    %   grid_num: Number of angle grid points for estimation (integer)
    %   is_plot: Flag to enable plotting (boolean)
    % Outputs: 
    %   Y_AoA: Updated struct with added field:
    %       - AoA_est: Estimated AoAs [L x 1] (degrees)
    % Example:
    %   [Y_AoA] = Func_MUSIC_AoA_Smooth(Y_AoA, sig_param, 1024, true)
    
    % Extract parameters from input structures
    Snap    = sig_param.Snap_AoA;
    Nr      = sig_param.Nr;
    K       = sig_param.K;
    Ls      = sig_param.Ls;
    Ld      = sig_param.Ld;
    W       = sig_param.W;
    AoA     = sig_param.AoA;
    b       = sig_param.b;
    
    Y       = Y_AoA.Y;
    L       = Ls + Ld;  % Total number of paths (static + dynamic)
    
    % Rearrange dimensions for processing: [Snap x Nr x K] -> [Nr x Snap x K]
    Y       = permute(Y, [2, 1, 3]);

    % Initialize correlation matrix
    R       = zeros(Nr, Nr);

    % Compute frequency smoothed correlation matrix averaged over all subcarriers
    % SEE: Equ. (14)
    for k = 1:K
        Y_temp = W * squeeze(Y(:, :, k));
        R      = R + Y_temp * Y_temp' / Snap;
    end
    
    % Normalize by number of subcarriers
    R = R / K;
    
    % Signal subspace (corresponding to L largest eigenvalues)
    [Us_AoA, eigvals] = eigs(R, L, 'largestabs');

    % Define angular grid for spectrum calculation
    angle_dist = linspace(-90, 90, grid_num);
    
    % Generate steering vectors for all angles in the grid
    ar = steervec(0.5 * (0:Nr-1), angle_dist);
    v  = ar;
    
    % Calculate MUSIC spectrum using noise subspace projection
    % SEE: Equ. (15)
    P_MUSIC_AoA = 1 ./ (Nr - sum(abs(v' * Us_AoA) .^ 2, 2));

    % Plot MUSIC spectrum and actual/estimated AoAs if requested
    if is_plot
        figure;
        semilogy(angle_dist, P_MUSIC_AoA, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
        hold on;
        
        % Plot static path (blue) and dynamic paths (red)
        stem(AoA(1), abs(b(1))/2, 'bo', 'MarkerSize', 10, 'LineWidth', 1.5);
        stem(AoA(2:end), abs(b(2:end))/2, 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
        
        Customed_Figure('Angle (degree)', 'P', 'MUSIC Spectrum for AoA Estimation');
    end

    % Find the peaks of P_MUSIC_AoA, which correspond to the estimated AoAs
    [pks, locs] = findpeaks(P_MUSIC_AoA);
    [pks, idx]  = sort(pks, 'descend');
    locs        = locs(idx);
    
    % Select the L highest peaks as the estimated AoAs
    AoA_est     = angle_dist(locs(1:L));
    
    % Display actual and estimated AoAs if plotting is enabled
    if is_plot
        disp("---------------------- AoA estimation ----------------------");
        disp(sprintf("The actual AoA is %f\n", sort(AoA, 'ascend')));
        disp(sprintf("The estimated AoA is %f\n", sort(AoA_est, 'ascend')));
        % disp("------------------------------------------------------------");
    end

    % Update output structure with estimated AoAs
    Y_AoA.AoA_est = AoA_est;
end