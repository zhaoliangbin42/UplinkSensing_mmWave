%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Function for calculating the array response using FFT
%                 Used to detect the Line-of-Sight (LoS) path and visualize
%                 array beam patterns in the angular domain
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rlt = Func_Array_Response(a, FFT_length, is_plot)
    % This function computes the array response in the angular domain using FFT
    %
    % Inputs:
    % a:           steering vector matrix (each column represents one steering vector)
    % FFT_length:  Length of the FFT (higher length means higher resolution)
    % is_plot:     Boolean flag to control visualization
    %
    % Output:
    % rlt:         Angular domain representation of the steering vectors
    
    % Normalize each steering vector column to have unit Frobenius norm
    for i = 1:size(a, 2)
        a(:, i) = a(:, i) / norm(a(:, i), 'fro');
    end
    
    % * NOTE: Alternative methods commented out for reference
    % a = a ./ max(abs(a));   % Normalize by maximum amplitude
    % a = a / norm(a);        % Normalize the entire matrix
    % a = [a(:); zeros(FFT_length - length(a), 1)]; % Zero padding
    % DFT_matrix = dftmtx(FFT_length);
    % rlt = DFT_matrix * a / FFT_length;
    
    % Perform FFT on each column of the steering vector matrix
    rlt = fft(a, FFT_length);
    
    % Shift zero frequency to center for proper visualization
    rlt = fftshift(rlt);

    % Create the angle codebook for the x-axis of the plot
    u_codebook = (-FFT_length/2 : FFT_length/2 - 1) / FFT_length;
    
    % Convert to angles in degrees using arcsine transformation
    % The factor of 2 scales the u_codebook range (-0.5,0.5) to (-1,1)
    angle_codebook = asin(u_codebook * 2) * 180 / pi;
    
    % Generate plot if requested
    if is_plot
        figure;
        
        % Plot each steering vector response
        for i = 1:size(a, 2)
            plot(angle_codebook, abs(rlt(:, i)), 'LineWidth', 2);
            
            % Find and mark the peak value
            [max_val, max_idx] = max(abs(rlt(:, i)));
            hold on;
            plot(angle_codebook(max_idx), max_val, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            text(angle_codebook(max_idx), max_val, num2str(max_val), 'FontSize', 16);
        end
        
        hold off;
        grid on;
        xlabel('Angle (degree)');
        ylabel('Amplitude');
        title('Array Response');

        % Apply custom figure formatting
        Customed_Figure('Angle (degree)', 'Amplitude', 'Array Response', [1 1 25 7]);
    end
end