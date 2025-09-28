clc; clear; close all;

%% 数据
% 地区名称
regions = {'地区A', '地区B', '地区C'};
% 地区A的数据
data_A = [
    47.22, 5.06, 0.10;
    47.45, 4.35, 0.15;
    47.52, 6.85, 0.12;
    47.86, 4.19, 0.17;
    47.31, 7.57, 0.18
];

% 地区B的数据
data_B = [
    54.33, 6.22, 0.12;
    56.17, 3.31, 0.15;
    48.40, 2.43, 0.22;
    52.62, 5.92, 1.12;
];

% 地区C的数据
data_C = [
    43.12, 10.33, 0.05;
    42.05, 9.67, 0.08;
    42.50, 9.62, 0.02;
    40.77, 9.68, 0.04
];

% 样本量
n_A = size(data_A, 1);
n_B = size(data_B, 1);
n_C = size(data_C, 1);
n_total = n_A + n_B + n_C;

% 成分名称
components = {'SiO_2', 'Fe_2O_3', 'K_2O'};

% 计算各组样本协方差矩阵
S_A = cov(data_A);
S_B = cov(data_B);
S_C = cov(data_C);

%% H0: Σ1 = Σ2 = Σ3 vs H1: 存在i≠j使得Σi ≠ Σj
fprintf('问题1：协方差矩阵相等性检验\n');
fprintf('H0: Σ1 = Σ2 = Σ3 vs H1: 存在i≠j使得Σi ≠ Σj\n');

% 基本参数
p = size(data_A, 2); % 变量个数 p=3
r = 3;               % 组数
n = n_total;         % 总样本数

% 各组样本协方差矩阵的行列式
det_S_A = det(S_A);
det_S_B = det(S_B);
det_S_C = det(S_C);

% 计算L矩阵
L_A = (n_A - 1) * S_A; % L_k = (n_k - 1) * S_k
L_B = (n_B - 1) * S_B;
L_C = (n_C - 1) * S_C;
L = L_A + L_B + L_C; 

% 合并协方差矩阵估计
S_pooled = L / (n - r);
det_S_pooled = det(S_pooled);


% 统计量
M = (n - r) * log(det(L)/(n-r)) - ...
    (n_A-1)*log(det(L_A)/(n_A-1)) - ...
    (n_B-1)*log(det(L_B)/(n_B-1)) - ...
    (n_C-1)*log(det(L_C)/(n_C-1));

fprintf('\n统计量计算：\n');
fprintf('M = %.6f\n', M);

% 计算参数
% 计算d1
sum_inv_ni_minus_1 = 1/(n_A-1) + 1/(n_B-1) + 1/(n_C-1);
inv_n_minus_r = 1/(n-r);
d1 = (2*p^2 + 3*p - 1)/(6*(p+1)*(r-1)) * (sum_inv_ni_minus_1 - inv_n_minus_r);

% 计算d2  
sum_inv_ni_minus_1_sq = 1/(n_A-1)^2 + 1/(n_B-1)^2 + 1/(n_C-1)^2;
inv_n_minus_r_sq = 1/(n-r)^2;
d2 = (p-1)*(p+2)/(6*(r-1)) * (sum_inv_ni_minus_1_sq - inv_n_minus_r_sq);

% 计算f1, f2, b
f1 = p*(p+1)*(r-1)/2;
f2 = (f1 + 2)/(d2 - d1^2);
b = f1 / (1 - d1 - f1/f2);

%% 统计检验
fprintf('\n假设检验：\n');

% 当样本量不是很大时，根据f2的值选择分布近似（查阅参考文献得知）
if f2 > 4
    % 使用F分布近似
    F_stat = M / b;
    df1 = f1;
    df2 = f2;
    
    fprintf('使用F分布近似：M ≈ b·F(f1,f2)\n');
    fprintf('F统计量 = M/b = %.6f\n', F_stat);
    fprintf('自由度：df1 = %d, df2 = %.4f\n', df1, df2);
    
    % F分布临界值
    F_critical = finv(0.95, df1, df2);
    fprintf('F临界值(α=0.05) = %.6f\n', F_critical);
    fprintf('p值 = %.6f\n', 1 - fcdf(F_stat, df1, df2));
    
    fprintf('\n决策：\n');
    if F_stat > F_critical
        fprintf('F = %.6f > %.6f，拒绝H0\n', F_stat, F_critical);
        fprintf('结论：在α=0.05水平下，三个地区的协方差矩阵不相等\n');
    else
        fprintf('F = %.6f ≤ %.6f，接受H0\n', F_stat, F_critical);
        fprintf('结论：在α=0.05水平下，三个地区的协方差矩阵相等\n');
    end

else
    % 使用卡方分布近似
    chi2_stat = M_alt * (1 - d1);
    df_chi2 = f1;
    
    fprintf('使用卡方分布近似：\n');
    fprintf('卡方统计量 = M(1-d1) = %.6f\n', chi2_stat);
    fprintf('自由度 = %.0f\n', df_chi2);
    
    chi2_critical = chi2inv(0.95, df_chi2);
    fprintf('卡方临界值(α=0.05) = %.6f\n', chi2_critical);
    fprintf('p值 = %.6f\n', 1 - chi2cdf(chi2_stat, df_chi2));
    
    fprintf('\n决策：\n');
    if chi2_stat > chi2_critical
        fprintf('χ² = %.6f > %.6f，拒绝H0\n', chi2_stat, chi2_critical);
        fprintf('结论：在α=0.05水平下，三个地区的协方差矩阵不相等\n');
    else
        fprintf('χ² = %.6f ≤ %.6f，接受H0\n', chi2_stat, chi2_critical);
        fprintf('结论：在α=0.05水平下，三个地区的协方差矩阵相等\n');
    end
end