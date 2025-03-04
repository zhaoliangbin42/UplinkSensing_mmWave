function Customed_Figure(customed_xlabel, customed_ylabel, customed_title, is_grid_on)

    if nargin < 5
        is_grid_on = true;
    end
    FontSize = 20;
    % 设置更高的输出分辨率
    % set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', customed_size);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);

    % 设置坐标轴标签和标题
    xlabel(customed_xlabel, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    ylabel(customed_ylabel, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    title(customed_title, 'FontName', 'Times New Roman', 'FontSize', FontSize);

    % 去掉外框
    % box off;

    % 开启网格
    if is_grid_on
        grid on;
    end

    % 优化图例和网格样式
    set(gca, 'GridLineStyle', '--', 'LineWidth', 1);
end