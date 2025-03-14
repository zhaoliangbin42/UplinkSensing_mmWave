%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Main execution script.
%
%    Tool versions: Matlab 2025a
%    Last update: 2025-03-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% Add required function paths
addpath("Funcs");
addpath("Gen_Data");

% Load global parameters
run('Params.m');

%% Configuration parameters
% Basic configuration
Ls            = 2;              % Number of static paths
Ld            = 3;              % Number of dynamic paths
L             = Ld + Ls;        % Total number of paths
EbN0_dB       = 10;             % Signal to noise ratio in dB
grid_num_AoA  = 1e5;            % Number of grids for AoA estimation
grid_num_DD   = 2048;           % Number of grids for Doppler-delay estimation
BF_type       = 'SINR';         % Beamforming type: 'SINR', 'NULL', 'BART', 'HYBD'
BF_rho        = 0.5;            % Weighting factor for hybrid beamforming
disp_pixel    = 1024;           % Number of pixels for display, disp_pixel x disp_pixel
rand('seed', 3);                % Set random seed for reproducibility

% Visualization control flags
is_plot_AoA_Est     = 1;              % Plot AoA estimation results
is_plot_LoS_Det     = 1;              % Plot Line-of-Sight path detection results
is_plot_BF_Response = 1;              % Plot combining vector response
is_plot_DD_Est      = 1;              % Plot Doppler and delay estimation results

% Pack parameters into structure for function calls
sig_param.Ld        = Ld;
sig_param.Ls        = Ls;
sig_param.L         = Ls + Ld;
sig_param.EbN0_dB   = EbN0_dB;

%% Generate random path parameters
% Generate distance, velocity, delay and Doppler for each path
% NOTE: First Ls paths are static (zero velocity), remaining Ld paths are dynamic
Distance_list = sort(rand(L, 1) * sig_param.max_distance, 'ascend'); % Distance [m], sorted in ascending order
Velocity_list = [zeros(Ls, 1); (rand(Ld, 1) - 1/2) * sig_param.max_velocity]; % Velocity [m/s]
Delay_list    = Distance_list / c;                  % Delay [s]
Doppler_list  = Velocity_list / c * fc;             % Doppler frequency [Hz]

% Generate Angles of Arrival (AoA) for all paths
% SEE: Method uses uniform spacing with random offset to avoid grid mismatch issues
AoA = (randperm(floor(2 * max_angle/min_interval), L) * min_interval - max_angle).';
AoA_offset = rand(1) * min_interval;
AoA = AoA + AoA_offset;
AoA = max(min(AoA, max_angle), -max_angle);         % Constrain AoA within allowed range
Ar = steervec(0.5 * [0:Nr - 1], AoA.');             % Generate steering vectors based on AoA

% Generate random complex path gains
% First path (assumed LoS) has unit magnitude with random phase
% Other paths have random complex gains with smaller magnitudes
b = [exp(1j * 2 * pi * rand(1)); 
     0.5 * (rand(Ls - 1, 1) + 1j * rand(Ls - 1, 1)); 
     0.5 * (rand(Ld, 1) + 1j * rand(Ld, 1))];

% Store all parameters in sig_param structure
sig_param.Delay_list   = Delay_list;
sig_param.Doppler_list = Doppler_list;
sig_param.AoA          = AoA;
sig_param.b            = b;

%% AES Stage
% Generate received signal for AoA estimation
Y_AoA = Gen_Y_AoA(sig_param); 

% Apply MUSIC algorithm for AoA estimation with frequency smoothing
% SEE: Section III-B
Y_AoA = Func_MUSIC_AoA_Smooth(Y_AoA, sig_param, grid_num_AoA, is_plot_AoA_Est);

% Detect Line-of-Sight path from estimated AoA results
% SEE: Section III-C
Y_AoA = Func_Detect_LoS(Y_AoA, sig_param, 1024, is_plot_LoS_Det);

%% DDE&CS Stage
% Step 1: BF vector optimization
% Optimize combining vector
% SEE: Section IV-A
w_opt = Func_Beamforming_Vec(Y_AoA, sig_param, BF_type, is_plot_BF_Response, BF_rho);

%% Step 2: Generate data for Doppler-delay estimation
% SEE: Section IV-B
% Generate signal for Doppler-delay estimation using optimized combining vector
Y_DD = Gen_Y_DD(sig_param, w_opt);

% Display ground truth Doppler and delay values
if is_plot_DD_Est
    disp('--------------- True Doppler and delay (unordered) ---------------');
    fprintf('Delay\t\t\t Doppler\n')
    for i = 1:L
        fprintf('%e s\t\t %e Hz\n', Delay_list(i), Doppler_list(i));
    end
end

%% Step 3: AoA-Based 2D-FFT-MUSIC (AB2FM) Algorithm
% SEE: Section IV-C, Algorithm 2
[delay_est, doppler_est] = Func_MUSIC_Doppler_Delay(Y_DD, sig_param, grid_num_DD, disp_pixel, is_plot_DD_Est, BF_type);