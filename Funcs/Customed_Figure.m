%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Custom MATLAB figure styling function with consistent 
%                 formatting
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Customed_Figure(customed_xlabel, customed_ylabel, customed_title, is_grid_on)
% CUSTOMED_FIGURE - Sets consistent styling for MATLAB figures
% 
% Applies a standardized style to the current figure with Times New Roman font,
% customized labels, and grid settings suitable for academic publications.
%
% Inputs:  
%   customed_xlabel: X-axis label text (string)
%   customed_ylabel: Y-axis label text (string)
%   customed_title:  Figure title text (string)
%   is_grid_on:      Flag to enable/disable grid (boolean, optional, default: true)
%
% Example:
%   Customed_Figure('Time (s)', 'Amplitude (dB)', 'Signal Response', true);

    % Default grid setting if not provided
    if nargin < 4
        is_grid_on = true;
    end
    
    % Define consistent font size for all text elements
    FontSize = 20;
    
    % Apply font settings to current axes
    set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);

    % Set axis labels and title with consistent formatting
    xlabel(customed_xlabel, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    ylabel(customed_ylabel, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    title(customed_title, 'FontName', 'Times New Roman', 'FontSize', FontSize);

    % NOTE: The following line is commented out but can be uncommented 
    % to remove the box around the plot
    % box off;

    % Apply grid settings based on input parameter
    if is_grid_on
        grid on;
        
        % Set grid style to dashed lines with specified width
        set(gca, 'GridLineStyle', '--', 'LineWidth', 1);
    end
end