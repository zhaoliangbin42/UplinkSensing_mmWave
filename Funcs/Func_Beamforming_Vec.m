%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Beamforming vector generation for mmWave ISAC systems
%                 with different optimization criteria (Bartlett beamformer,
%                 SINR-maximizing beamformer, null-space approach, and hybrid approach).
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_opt = Func_Beamforming_Vec(Y_AoA, sig_param, BF_type, is_plot, rho)
% FUNC_BEAMFORMING_VEC - Generate optimal beamforming vectors for various criteria
%
% Inputs:
%   Y_AoA:     Received signal structure containing:
%              - Y: Received signal matrix (Nr x Snap_AoA x K)
%   sig_param: Signal parameters structure containing:
%              - Snap_AoA: Number of snapshots for AoA estimation
%              - Nr: Number of receive antennas
%              - K: Number of subcarriers
%              - W: Combining matrix of the angle codebook
%              - L: Number of paths
%              - AoA: Actual Angles of Arrival (degrees)
%              - b: Path gains for verification
%   BF_type:   Beamforming type (string): 'SINR', 'NULL', 'BART', or 'HYBD'
%              Options:
%              'SINR': SINR-maximizing beamformer
%              'NULL': Null-space approach
%              'BART': Bartlett beamformer
%              'HYBD': Hybrid beamformer, combining Bartlett and NULL with weighting factor 'rho'
%   is_plot:   Flag to enable/disable plotting results (boolean)
%   rho:       Weighting factor for hybrid beamforming (0-1), 
%              the lower the value, the higher the interference suppression
%
% Outputs:
%   w_opt:     Optimal beamforming vector (Nr x 1)
%
% Example:
%   w_opt = Func_Beamforming_Vec(Y_AoA, sig_param, 'SINR', true)

    if nargin < 5
        rho = 0.05;  % Default weighting factor for hybrid beamforming
    end

    % Validate beamforming type
    BF_type = lower(BF_type);
    valid_bf_types = ["sinr", "null", "bart", 'hybd'];
    if ~ismember(BF_type, valid_bf_types)
        error("Invalid beamforming type, please choose from 'SINR', 'NULL', 'BART', 'HYBD'");
    end

    % Extract parameters from input structures
    Y        = Y_AoA.Y;
    Snap_AoA = sig_param.Snap_AoA;
    Nr       = sig_param.Nr;
    AoA      = sig_param.AoA;
    W        = sig_param.W;
    L        = sig_param.L;

    % Generate steering vectors for all paths
    % NOTE: steervec expects array element spacing in wavelengths (0.5 = half-wavelength)
    Ar       = steervec(0.5 * [0:Nr-1], AoA.');
    Ar_LoS   = Ar(:, 1);          % LoS path steering vector
    Ar_NLoS  = Ar(:, 2:end);      % NLoS paths steering vectors
    
    % Process received signal
    % NOTE: Only using first subcarrier for beamforming optimization
    Y = squeeze(Y(:, :, 1)).';    % Dimensions: Snap_AoA x Nr
    H = W * Y;                    % Channel matrix in beamspace
    
    % Select beamforming method based on input parameter
    switch BF_type
        case 'sinr'
            % SINR-maximizing beamformer - optimizes signal-to-interference-plus-noise ratio
            w_opt = Func_Opt_SINR(H, Ar_LoS, Ar_NLoS);
        case 'null'
            % Null-space approach
            w_opt = Func_Opt_NULL(H, Ar_LoS, Ar_NLoS);
        case 'bart'
            % Bartlett beamformer
            w_opt = Ar_LoS;
        case 'hybd'
            % Hybrid beamformer - combines Bartlett and NULL with weighting factor
            w_opt = Func_Opt_Hybrid(H, Ar_LoS, Ar_NLoS, rho);
    end

    % Visualization if requested
    if is_plot
        % Plot array response pattern and mark actual AoA locations
        Func_Array_Response(w_opt, 2048, 1);
        hold on;
        
        % Highlight LoS path with blue marker
        stem(AoA(1), abs(sig_param.b(1)), 'bo', 'LineWidth', 3);
        
        % Mark NLoS paths with red markers
        stem(AoA(2:end), abs(sig_param.b(2:end)), 'ro', 'LineWidth', 1.5);
        
        title(['Beamforming Pattern: ' upper(BF_type)]);
    end

end