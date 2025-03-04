function Customed_Figure(customed_xlabel, customed_ylabel, customed_title, is_grid_on)

    if nargin < 5
        is_grid_on = true;
    end
    FontSize = 20;
    % ���ø��ߵ�����ֱ���
    % set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', customed_size);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize);

    % �����������ǩ�ͱ���
    xlabel(customed_xlabel, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    ylabel(customed_ylabel, 'FontName', 'Times New Roman', 'FontSize', FontSize);
    title(customed_title, 'FontName', 'Times New Roman', 'FontSize', FontSize);

    % ȥ�����
    % box off;

    % ��������
    if is_grid_on
        grid on;
    end

    % �Ż�ͼ����������ʽ
    set(gca, 'GridLineStyle', '--', 'LineWidth', 1);
end