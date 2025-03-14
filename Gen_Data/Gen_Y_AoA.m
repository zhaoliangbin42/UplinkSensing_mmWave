%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Generate received signal data for Angle of Arrival (AoA) estimation
%                 in mmWave Uplink Sensing systems. Models multi-path propagation with
%                 both static and dynamic paths, accounting for Doppler shifts, delays, 
%                 timing offsets, and carrier frequency offsets.
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y_AoA = Gen_Y_AoA(sig_param)
% GEN_Y_AOA - Generate received signal data for AoA estimation 
% 
% Inputs:  
%   sig_param: Structure containing system and signal parameters
%     - Nr: Number of receive antennas
%     - K: Number of subcarriers
%     - Ts: Sample time [s]
%     - Tm: Frame duration [s]
%     - fc: Carrier frequency [Hz]
%     - D_f: Subcarrier spacing [Hz]
%     - W: Beamforming codebook matrix [Nr x codebook_size]
%     - Snap_AoA: Number of AoA snapshots
%     - Ld: Number of dynamic paths
%     - Ls: Number of static paths
%     - EbN0_dB: Signal-to-noise ratio [dB]
%     - AoA: Angle of arrival list [rad]
%     - b: Path gain coefficients [complex]
%     - Delay_list: Path delay list [s]
%     - Doppler_list: Doppler frequency list [Hz]
%
% Outputs: 
%   Y_AoA: Structure containing generated received signals
%     - Y: Received signal for each snapshot [Snap_AoA x Nr x K]
%
% Example:
%   Y_AoA = Gen_Y_AoA(sig_param)
%
% SEE: Multi-snapshot implementation of Equ. (5)

    % Extract parameters from input structure
    Nr         = sig_param.Nr;          % Number of receive antennas
    K          = sig_param.K;           % Number of subcarriers
    Ts         = sig_param.Ts;          % Sample time [s]
    Tm         = sig_param.Tm;          % Frame duration [s]
    fc         = sig_param.fc;          % Carrier frequency [Hz]
    D_f        = sig_param.D_f;         % Subcarrier spacing [Hz]
    Snap_AoA   = sig_param.Snap_AoA;    % Number of snapshots for AoA estimation
    Ld         = sig_param.Ld;          % Number of dynamic paths
    Ls         = sig_param.Ls;          % Number of static paths
    L          = Ld + Ls;               % Total number of paths
    EbN0_dB    = sig_param.EbN0_dB;     % Signal-to-noise ratio [dB]
    AoA        = sig_param.AoA;         % Angle of arrival list [rad]
    b          = sig_param.b;           % Path gain coefficients
    W          = sig_param.W;           % Beamforming codebook
    Delay_list = sig_param.Delay_list;  % Path delay list [s]
    Doppler_list = sig_param.Doppler_list; % Doppler frequency list [Hz]

    % Initialize output structure
    Y_AoA = struct();

    % Generate steering vectors for each AoA
    % SEE: Equ. (3)
    Ar = steervec(0.5 * (0:Nr-1), AoA.');

    % Pre-allocate signal arrays
    Y            = zeros(Snap_AoA, Nr, K);  % Received signal array
    
    % Subcarrier frequency vector
    subcarrier_freqs = fc + (1:K) * D_f;
    
    % Calculate noise power from SNR
    noise_power = 10^(-EbN0_dB/20) / sqrt(2);
    
    % Generate snapshots
    for i = 1:Snap_AoA 
        % Generate random time and carrier frequency offsets
        TO  = rand(1) * 2 * pi;  % Time offset [0, 2¦Ð]
        CFO = rand(1) * 2 * pi;  % Carrier frequency offset [0, 2¦Ð]
        
        % Process each antenna element
        for j = 1:Nr
            % Calculate time-dependent Doppler shift
            % NOTE: Doppler varies with snapshot and antenna index
            time_index = (i-1) * Tm + j;
            Dopp = Doppler_list * Ts * time_index;
            
            % Generate complex Gaussian noise
            n = (randn(Nr, K) + 1j * randn(Nr, K)) * noise_power;
            
            % Calculate time and carrier frequency offset factor
            TO_CFO = exp(-1j * 2 * pi * TO * subcarrier_freqs) * ...
                     exp(1j * 2 * pi * CFO * Ts * time_index);
            
            % Compute channel response incorporating:
            % - Time/frequency offsets 
            % - Path gains with Doppler shifts
            % - Frequency-dependent phase shifts due to delay
            % - AWGN noise
            % SEE: Equ. (6)(7)(8)
            H = TO_CFO .* (Ar * (b .* exp(1j * 2 * pi * Dopp) .* ...
                exp(-1j * 2 * pi * Delay_list .* subcarrier_freqs))) + n;
            
            % Apply beamforming weight vector
            % SEE: Equ. (5)
            Y(i, j, :) = W(:, j)' * H;
        end
    end

    % Store outputs
    Y_AoA.Y      = Y;
    
end