%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Benko
%    Email: lb_zhao_bit_ee@163.com
%
%    Description: Optimal combining vector computation for nulling NLoS paths
%                 and maximizing SNR in the remaining signal subspace.
%
%    Tool versions: Matlab 2025a
%    Affiliation: Beijing Institute of Technology
%    Last update: 2025-03-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w_opt = Func_Opt_NULL(H, Ar_LoS, Ar_NLoS)
    % FUNC_OPT_NULL - Computes optimal combining vector by nulling NLoS paths
    % Inputs:
    %   H: Channel matrix (complex matrix, size [M x Snapshot])
    %   Ar_LoS: Line-of-Sight steering vector (complex vector, size [M x 1])
    %   Ar_NLoS: Non-Line-of-Sight steering vectors (complex matrix, size [M x (L-1)])
    % Outputs:
    %   w_opt: Optimal combining vector (complex vector, size [M x 1])


    % Derive parameters from input matrices
    Snapshot = size(H, 2);    % Number of AoA snapshots
    L        = size(Ar_NLoS, 2) + 1; % Total number of paths (LoS + NLoS)
    Ar       = [Ar_LoS, Ar_NLoS];    % Combined steering vectors

    % Compute the nullspace of the NLoS paths
    % NOTE: This creates a basis for the subspace orthogonal to all NLoS paths
    W_NLoS_null = null(Ar_NLoS.'); 

    % Project the channel matrix onto the nullspace of NLoS paths
    % This effectively removes the NLoS components from the channel
    NH = W_NLoS_null.' * H; 

    % Find the optimal combining vector using eigendecomposition
    % This maximizes the received signal power in the nullspace of NLoS paths
    % SEE: Equ. (23) and around
    [eigvec, eigval] = eig(NH * NH');
    eigval = diag(eigval);  % Extract diagonal elements from eigenvalue matrix
    [~, idx] = sort(eigval, 'descend');  % Sort eigenvalues in descending order

    % Select the eigenvector corresponding to the largest eigenvalue
    % and project it back to the original signal space
    w_opt = eigvec(:, idx(1)).' * W_NLoS_null';
    w_opt = w_opt(:);  % Ensure column vector output
end