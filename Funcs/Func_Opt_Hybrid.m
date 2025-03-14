%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Hybrid combining vector computation that combines 
%                 Bartlett beamformer and null-space approach with a tunable parameter.
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_opt = Func_Opt_Hybrid(H, Ar_LoS, Ar_NLoS, rho)
    % FUNC_OPT_HYBRID - Computes hybrid combining vector using weighted combination 
    %                   of Bartlett beamformer and null-space approach
    % Inputs:  
    %   H: Channel matrix (complex matrix, size [M x Snapshot])
    %   Ar_LoS: Line-of-Sight steering vector (complex vector, size [M x 1])
    %   Ar_NLoS: Non-Line-of-Sight steering vectors (complex matrix, size [M x (L-1)])
    %   rho: Weighting factor between Bartlett beamformer and null-space approach (scalar, 0 ¡Ü rho ¡Ü 1)
    %        the lower the value, the higher the interference suppression
    % Outputs: 
    %   w_opt: Optimal hybrid combining vector (complex vector, size [M x 1])
    % Example:
    %   w_opt = Func_Opt_Hybrid(H, Ar_LoS, Ar_NLoS, 0.1)


    % Derive parameters from input matrices
    Snapshot = size(H, 2);    % Number of AoA snapshots
    L        = size(Ar_NLoS, 2) + 1; % Total number of paths (LoS + NLoS)
    Ar       = [Ar_LoS, Ar_NLoS];    % Combined steering vectors

    % Compute combining vector based on null-space approach
    w_null = Func_Opt_NULL(H, Ar_LoS, Ar_NLoS);

    % Calculate phase difference between NS vector and Bartlett beamformer
    % NOTE: This phase alignment ensures constructive combination of vectors
    phi = angle(w_null' * Ar_LoS);

    % Create hybrid combining vector as weighted sum of Bartlett beamformer 
    % and phase-adjusted NS vector
    % SEE: Equ. (24)
    w_opt = sqrt(rho) .* Ar_LoS + sqrt(1 - rho) .* w_null * exp(1j * phi);
    
    % Normalize the combining vector to unit norm
    % NOTE: This maintains constant transmit power constraint
    w_opt = w_opt ./ vecnorm(w_opt);
end