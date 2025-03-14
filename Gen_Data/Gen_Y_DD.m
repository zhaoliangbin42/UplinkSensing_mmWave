%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Generate data for AB2FM algorithm simulation
%                 with time offset and carrier frequency offset compensation
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zeta_k = Gen_Y_DD(sig_param, w_ref)
    % GEN_Y_DD - Generates signals in delay-Doppler domain with TO/CFO compensation
    % Inputs:  
    %   sig_param: Structure with signal parameters
    %     - Snap_DD: Number of snapshots
    %     - Nr: Number of receive antennas
    %     - K: Number of subcarriers
    %     - M: Number of OFDM symbols per subcarrier
    %     - Ts: Sampling time (seconds)
    %     - Tm: OFDM symbol duration (seconds)
    %     - fc: Carrier frequency (Hz)
    %     - D_f: Subcarrier spacing (Hz)
    %     - EbN0_dB: Energy per bit to noise ratio (dB)
    %     - Doppler_list: List of Doppler shifts for each path (Hz)
    %     - Delay_list: List of delays for each path (seconds)
    %     - AoA: Angles of arrival (radians)
    %     - b: Path gains (complex)
    %     - W: Combining vector for reception
    %     - L: Number of multipath components
    %   w_ref: Reference weight vector for TO/CFO compensation
    % Outputs: 
    %   zeta_k: Compensated received signal [Nr¡ÁK¡ÁSnap_DD]

    % Extract parameters from input structure
    Snap_DD      = sig_param.Snap_DD;
    Nr           = sig_param.Nr;
    K            = sig_param.K;
    M            = sig_param.M;
    Ts           = sig_param.Ts;
    Tm           = sig_param.Tm;
    fc           = sig_param.fc;
    D_f          = sig_param.D_f;
    EbN0_dB      = sig_param.EbN0_dB;
    Doppler_list = sig_param.Doppler_list;
    Delay_list   = sig_param.Delay_list;
    AoA          = sig_param.AoA;
    b            = sig_param.b;
    W            = sig_param.W;
    L            = sig_param.L;

    % Generate steering vectors for all antennas and angles
    Ar = steervec(0.5 * (0:Nr-1), AoA.');  % Steering vector matrix [Nr¡ÁL]
    
    % Pre-allocate received signal matrix
    % NOTE: zeta_k represents the received signal after TO/CFO compensation
    zeta_k = zeros(Nr, K, Snap_DD);  
    
    % Process each snapshot
    for i = 1:Snap_DD
        for j = 1:Nr
            % Generate random time and frequency offsets
            TO  = rand(1) * 2 * pi;  % Random time offset
            CFO = rand(1) * 2 * pi;  % Random carrier frequency offset

            % Generate random BPSK symbols {-1,+1}
            s = randi([0 1], 1, M*K) * 2 - 1;
            
            % Calculate Doppler effect for each path over time
            time_idx = (i * Nr + j) * Tm + (1:M);
            Doppler = repmat(Doppler_list * Ts * time_idx, 1, K);  % Doppler phase for each path
            
            % Calculate delay effect for each path over frequency
            freq_idx = fc + (1:K) * D_f;
            Delay = repelem(Delay_list * freq_idx, 1, M);  % Delay phase for each path
            
            % Generate complex AWGN noise for two branches (reference and data)
            noise_factor = 10^(-EbN0_dB/20);
            n1 = (randn(Nr, M*K) + 1j * randn(Nr, M*K)) / sqrt(2) * noise_factor;
            n2 = (randn(Nr, M*K) + 1j * randn(Nr, M*K)) / sqrt(2) * noise_factor;

            % Generate time and frequency offset effects
            % NOTE: Combined TO and CFO for all subcarriers and symbols
            TO_phase = -1j * 2 * pi * TO * freq_idx;
            CFO_phase = 1j * 2 * pi * CFO * Ts * time_idx.';
            TO_CFO_base = exp(TO_phase) .* exp(CFO_phase);
            TO_CFO = repmat(TO_CFO_base(:).', L, 1);

            % Generate channel matrix with Doppler and delay effects
            channel_phases = TO_CFO .* (b .* exp(1j * 2 * pi * Doppler) .* exp(-1j * 2 * pi * Delay));
            H = Ar * (channel_phases .* s);  % Apply steering vectors and symbols
            
            % Add noise to create two noisy channel observations
            H1 = H + n1;  % For data path
            H2 = H + n2;  % For reference path
            
            % Apply combining vectors
            Y_s = W(:, j)' * H1;    % Data signal, SEE Equ. (30)
            Y_c = w_ref' * H2;    % Reference signal, SEE Equ. (29)
            
            % Compensate for TO/CFO
            % SEE: Equ. (34)
            compensated = Y_s ./ Y_c;
            zeta_k(j, :, i) = reshape(mean(reshape(compensated.', M, []), 1), K, []).';
        end
    end

end