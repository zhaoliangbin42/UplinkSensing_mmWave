%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Author: Liangbin
%    Email: zhaoliangbin@bit.edu.cn
%
%    Description: Function for comparing two sets of numbers with tolerance
%
%    Tool versions: Matlab 2025a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function is_equal_num = Func_Compare_Sets(set1, set2, tolerance)
    % FUNC_COMPARE_SETS - Compares elements between two sets within a specified tolerance
    % Inputs:  
    %   set1: First set of numeric values (array)
    %   set2: Second set of numeric values (array)
    %   tolerance: Maximum acceptable difference between matching elements (scalar)
    % Outputs: 
    %   is_equal_num: Number of elements in set1 that match with any element in set2 (integer)
    % Example:
    %   matches = Func_Compare_Sets([10.1, 20.2], [10.0, 20.3], 0.5)

    % Validate input arguments
    if length(set1) ~= length(set2)
        error('The two sets must be of the same length.');
    end
    
    % Ensure column vectors for consistent processing
    set1 = set1(:);
    set2 = set2(:);
    
    % Pre-allocate result array for matched elements
    compare_list = zeros(length(set1), 1);  % Pre-allocation
    
    % For each element in set1, check if there's a matching element in set2 within tolerance
    for i = 1:length(set1)
        % Calculate absolute differences between current element and all elements in set2
        diff_values  = abs(set1(i) - set2);
        
        % Find the minimum difference
        min_diff     = min(diff_values);
        
        % Check if the minimum difference is within tolerance
        is_matched   = (min_diff < tolerance);
        
        % Record the match result
        compare_list(i) = is_matched;
    end
    
    % NOTE: Alternative vectorized implementation (commented out)
    % This approach assumes elements at the same position should match
    % set1_sort = sort(set1, 'ascend');
    % set2_sort = sort(set2, 'ascend');
    % compare_list = abs(set1_sort - set2_sort) < tolerance;
    
    % Count the total number of matched elements
    is_equal_num = sum(compare_list);
end