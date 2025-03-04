%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Benko
%    Email: lb_zhao_bit_ee@163.com
%
%    Description: Optimal beamformer design based on SINR maximization
%                 for mmWave sensing with Line-of-Sight and Non-Line-of-Sight paths
%
%    Tool versions: Matlab 2025a
%    Affiliation: Beijing Institute of Technology
%    Last update: 2025-03-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_opt = Func_Opt_SINR(H, Ar_LoS, Ar_NLoS)
    % FUNC_OPT_SINR - Calculate optimal beamforming vector to maximize SINR
    % 
    % Inputs:
    %   H:       Channel measurements matrix (Nr x Snapshot)
    %   Ar_LoS:  Steering vector for Line-of-Sight path (Nr x 1)
    %   Ar_NLoS: Steering vectors for Non-Line-of-Sight paths (Nr x L-1)
    %
    % Outputs:
    %   w_opt:   Optimal beamforming vector that maximizes SINR (Nr x 1)

    % Derive parameters from input matrices
    Snapshot = size(H, 2);    % Number of AoA snapshots
    L        = size(Ar_NLoS, 2) + 1; % Total number of paths (LoS + NLoS)
    Ar       = [Ar_LoS, Ar_NLoS];    % Combined steering vectors

    % Noise power estimation using eigenvalue decomposition
    % NOTE: Assumes noise subspace corresponds to smallest eigenvalues
    [eigvec, eigval] = eig(H * H');
    eigval = diag(eigval);
    [~, idx] = sort(eigval, 'descend');
    N_power = mean(eigval(idx(L+1:end)));  % Average of noise subspace eigenvalues

    % Least square estimation of path gains
    % SEE: Equ. (26)
    p_A    = pinv(Ar);      % Pseudo-inverse of steering vector matrix
    b      = p_A * H;       % Estimated path gains (includes noise)
    b0     = b(1, :);       % LoS path gain
    bNLoS  = b(2:end, :);   % NLoS path gains

    % Construct signal covariance matrices
    R0 = Ar_LoS * b0 * b0' * Ar_LoS' / Snapshot;    % Covariance matrix of LoS path
    
    % NOTE: Create NLoS + noise covariance matrix
    RN = Ar_NLoS * bNLoS * bNLoS' * Ar_NLoS' / Snapshot + eye(size(Ar_NLoS, 1)) * N_power;

    % Find optimal beamforming vector using generalized eigenvalue problem
    % The optimal w maximizes w'*R0*w / w'*RN*w
    % SEE: Equ. (27)
    RR = pinv(RN) * R0;   % SINR optimization matrix
    [eigvec, eigval] = eig(RR);
    eigval = diag(eigval);
    [~, idx] = sort(eigval, 'descend');
    w_opt = eigvec(:, idx(1));  % Eigenvector corresponding to largest eigenvalue
    w_opt = w_opt(:);           % Ensure column vector format

end