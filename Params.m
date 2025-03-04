% User defined parameters
Nr       = 24;     % Number of antennas
Snap_AoA = 50;     % Number of snapshots for Angle of Arrival estimation
Snap_DD  = 10;     % Number of snapshots for Distance-Doppler estimation
K        = 32;     % Number of subcarriers
M        = 32;     % Length of each frame (symbols per frame)
Tm       = 1024;   % Interval between two frames (in samples)

% System configurations
fc       = 26e9;   % Carrier frequency, 26 GHz (mmWave band)
D_f      = 100e6 / K;  % Subcarrier spacing (Hz)
Ts       = 1 / D_f;    % Sampling interval (s), equal to symbol duration
c        = 3e8;    % Speed of light (m/s)

% Calculate sensing limitations based on system parameters
max_distance  = c / D_f;               % Maximum unambiguous distance (m)
max_velocity  = c / (Ts * Tm * fc);    % Maximum unambiguous velocity (m/s)

% Beamspace codebook generation
u_codebook = (0:Nr - 1) / Nr;    % Normalized spatial frequencies
% * NOTE: Wrap spatial frequencies to [-0.5, 0.5] range for beamforming
u_codebook(u_codebook > 0.5) = u_codebook(u_codebook > 0.5) - 1;
angle_codebook = asin(u_codebook * 2) * 180 / pi;  % Convert to physical angles (degrees)
W = steervec(0.5 * [0:Nr - 1], angle_codebook);    % Combining matrix of the angle codebook
W = W / norm(W);   % Normalize beamforming matrix

% Create parameter structure for easy passing to functions
sig_param = struct();
sig_param.Nr           = Nr;           % Number of receive antennas
sig_param.K            = K;            % Number of subcarriers
sig_param.M            = M;            % Length of each frame
sig_param.Tm           = Tm;           % Interval between two frames
sig_param.fc           = fc;           % Carrier frequency
sig_param.D_f          = D_f;          % Subcarrier spacing
sig_param.Ts           = Ts;           % Sampling interval, symbol rate
sig_param.c            = c;            % Speed of light
sig_param.W            = W;            % Combining matrix of the angle codebook
sig_param.Snap_AoA     = Snap_AoA;     % Number of snapshots for AoA estimation
sig_param.Snap_DD      = Snap_DD;      % Number of snapshots for Doppler and delay estimation
sig_param.max_distance = max_distance; % Maximum unambiguous distance
sig_param.max_velocity = max_velocity; % Maximum unambiguous velocity

% Parameters for AoA generation
max_angle    = 80;    % Maximum angle in degrees (range: -max_angle to +max_angle)
min_interval = 5;     % Minimum angular separation between targets (degrees)