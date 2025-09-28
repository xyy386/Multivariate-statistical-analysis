clc; clear; close all;

%% 数据
% 国家名称
countries = {'阿根廷', '澳大利亚', '奥地利', '比利时', '百慕大', ...
            '巴西', '缅甸', '加拿大', '智利', '中国', '哥伦比亚'};

% 原始数据矩阵 (11个国家 × 8个项目)
% 列顺序：100米(秒), 200米(秒), 400米(秒), 800米(分), 1500米(分), 5000米(分), 10000米(分), 马拉松(分)
raw_data = [
    10.39, 20.81, 46.84, 1.81, 3.7,  14.04, 29.36, 137.72;
    10.31, 20.06, 44.84, 1.74, 3.57, 13.28, 27.66, 128.3;
    10.44, 20.81, 46.82, 1.79, 3.6,  13.26, 27.72, 135.9;
    10.34, 20.68, 45.04, 1.73, 3.6,  13.22, 27.45, 129.95;
    10.28, 20.58, 45.91, 1.8,  3.75, 14.68, 30.55, 146.62;
    10.22, 20.43, 45.21, 1.73, 3.66, 13.62, 28.62, 133.13;
    10.64, 21.52, 48.3,  1.8,  3.85, 14.45, 30.28, 139.95;
    10.17, 20.22, 45.68, 1.76, 3.63, 13.55, 28.09, 130.15;
    10.34, 20.8,  46.2,  1.79, 3.71, 13.61, 29.3,  134.03;
    10.51, 21.04, 47.3,  1.81, 3.73, 13.9,  29.13, 133.53;
    10.43, 21.05, 46.1,  1.82, 3.74, 13.49, 27.88, 131.35
];

data = raw_data;
% 将分钟转换为秒（第4-8列）
data(:, 4:8) = data(:, 4:8) * 60;

% 计算样本均值向量
mu = mean(data);
fprintf('各项目平均成绩（秒）：\n');
events = {'100米', '200米', '400米', '800米', '1500米', '5000米', '10000米', '马拉松'};
for i = 1:length(events)
    if i <= 3
        fprintf('%s: %.2f秒\n', events{i}, mu(i));
    else
        fprintf('%s: %.2f秒 (%.2f分钟)\n', events{i}, mu(i), mu(i)/60);
    end
end
fprintf('\n');

% 计算协方差矩阵
S = cov(data);

%% 计算马氏距离
n = size(data, 1);
mahal_distances = zeros(n, 1);

for i = 1:n
    % 计算每个国家与均值的差向量
    diff_vec = data(i, :) - mu;
    
    % 计算马氏距离的平方：d²ₘ = (X - μ)' Σ⁻¹ (X - μ)
    mahal_distances(i) = diff_vec * inv(S) * diff_vec';
    
    fprintf('%s: d²ₘ = %.6f\n', countries{i}, mahal_distances(i));
end

%% 排序分析
fprintf('\n排序:\n');

% 按马氏距离排序
[sorted_dist_sq, sort_idx] = sort(mahal_distances);
fprintf('\n按马氏距离(d²ₘ)排序（从小到大）：\n');
fprintf('排名\t国家\td²ₘ\n');
for i = 1:n
    idx = sort_idx(i);
    fprintf('\n%d\t%s\t%.6f\t%.6f\n', i, countries{idx}, ...
            mahal_distances(idx));
end

%% 可视化分析
figure('Position', [100, 100, 1200, 400]);

% 图1：马氏距离柱状图
subplot(1, 3, 1);
bar(mahal_distances(sort_idx));
set(gca, 'XTickLabel', countries(sort_idx), 'XTickLabelRotation', 45);
title('各国马氏距离排序');
ylabel('马氏距离 d²ₘ');
grid on;

% 图2：原始数据热力图
subplot(1, 3, 2);
imagesc(data);
colorbar;
set(gca, 'XTick', 1:8, 'XTickLabel', events, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:11, 'YTickLabel', countries);
title('原始成绩数据热力图');

% 图3：标准化后的数据
subplot(1, 3, 3);
standardized_data = zscore(data);
imagesc(standardized_data);
colorbar;
set(gca, 'XTick', 1:8, 'XTickLabel', events, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:11, 'YTickLabel', countries);
title('标准化后数据热力图');

sgtitle('田径成绩马氏距离分析', 'FontSize', 16);
