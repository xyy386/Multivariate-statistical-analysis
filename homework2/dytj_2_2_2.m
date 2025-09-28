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

% 样本量
n_A = size(data_A, 1);
n_B = size(data_B, 1);
n_total = n_A + n_B;
% 成分名称
components = {'SiO_2', 'Fe_2O_3', 'K_2O'};

% 基本参数
p = size(data_A, 2); % 变量个数 p=3
r = 3;               % 组数
n = n_total;         % 总样本数

%% 问题2：两个正态总体均值相等性检验（地区A和B）
fprintf('问题2：地区A和B的均值相等性检验\n');

% 计算样本均值
X_bar_A = mean(data_A)'; % 列向量
X_bar_B = mean(data_B)'; % 列向量


fprintf('地区A样本均值：\n');
for i = 1:length(components)
    fprintf('  %s: %.4f\n', components{i}, X_bar_A(i));
end
fprintf('地区B样本均值：\n');
for i = 1:length(components)
    fprintf('  %s: %.4f\n', components{i}, X_bar_B(i));
end

% 计算A和B的协方差矩阵
S_A = cov(data_A);
S_B = cov(data_B);


%% Case1：假设协方差矩阵相等
fprintf('=================================')
fprintf('\nCase1：假设协方差矩阵相等的检验\n');
fprintf('H0: μ1 = μ2 vs H1: μ1 ≠ μ2\n');

% 合并协方差矩阵估计
S_pooled_AB = ((n_A - 1) * S_A + (n_B - 1) * S_B) / (n_A + n_B - 2);


% T²统计量计算（等方差情况）
X_diff = X_bar_A - X_bar_B;
T2_equal = (n_A * n_B) / (n_A + n_B) * X_diff' * inv(S_pooled_AB) * X_diff;

% F统计量转换
F_equal = ((n_A + n_B - p - 1) / (p * (n_A + n_B - 2))) * T2_equal;
df1_equal = p;
df2_equal = n_A + n_B - p - 1;

fprintf('\n统计量计算：\n');
fprintf('T² = %.6f\n', T2_equal);
fprintf('F = %.6f\n', F_equal);
fprintf('自由度：df1 = %d, df2 = %d\n', df1_equal, df2_equal);

% 临界值和p值
F_critical_equal = finv(0.95, df1_equal, df2_equal);
p_value_equal = 1 - fcdf(F_equal, df1_equal, df2_equal);

fprintf('F临界值(α=0.05) = %.6f\n', F_critical_equal);
fprintf('p值 = %.6f\n', p_value_equal);

fprintf('\n决策：\n');
if F_equal > F_critical_equal
    fprintf('F = %.6f > %.6f，拒绝H0\n', F_equal, F_critical_equal);
    fprintf('结论：在α=0.05水平下，地区A和B的均值向量显著不同\n');
else
    fprintf('F = %.6f ≤ %.6f，接受H0\n', F_equal, F_critical_equal);
    fprintf('结论：在α=0.05水平下，地区A和B的均值向量无显著差异\n');
end

%% Case2：协方差矩阵不相等的检验
fprintf('=================================')
fprintf('\nCase2：协方差矩阵不相等的检验\n');
fprintf('H0: μ1 = μ2 vs H1: μ1 ≠ μ2\n');

% 计算Lx和Ly矩阵
L_x = (n_A - 1) * S_A;
L_y = (n_B - 1) * S_B;

% 计算S*矩阵：S* = L_x/(n1(n1-1)) + L_y/(n2(n2-1))
S_star = L_x/(n_A*(n_A-1)) + L_y/(n_B*(n_B-1));

% 计算T²统计量
T2_BF = X_diff' * inv(S_star) * X_diff;

% 计算f^(-1)系数
% 第一项：(n1³-n1²)^(-1) * [(X̄-Ȳ)'S*^(-1)(L_x/(n1-1))S*^(-1)(X̄-Ȳ)]² * T^(-4)
term1_coeff = (n_A^3 - n_A^2)^(-1);
L_x_normalized = L_x / (n_A - 1);
quad_form1 = X_diff' * inv(S_star) * L_x_normalized * inv(S_star) * X_diff;
term1 = term1_coeff * (quad_form1^2) * (T2_BF^(-4));

% 第二项：(n2³-n2²)^(-1) * [(X̄-Ȳ)'S*^(-1)(L_y/(n2-1))S*^(-1)(X̄-Ȳ)]² * T^(-4)
term2_coeff = (n_B^3 - n_B^2)^(-1);
L_y_normalized = L_y / (n_B - 1);
quad_form2 = X_diff' * inv(S_star) * L_y_normalized * inv(S_star) * X_diff;
term2 = term2_coeff * (quad_form2^2) * (T2_BF^(-4));

f_inv = term1 + term2;
f_param = 1 / f_inv;

% 当H0成立时，((f-p+1)/fp)T² ~ F_{p,f-p+1}
if f_param > p
    F_BF = ((f_param - p + 1) / (f_param * p)) * T2_BF;
    df1_BF = p;
    df2_BF = f_param - p + 1;
    
    fprintf('\n统计量计算：\n');

    fprintf('T² = %.6f\n', T2_BF);
    fprintf('F = ((f-p+1)/(fp))T² = %.6f\n', F_BF);
    fprintf('自由度：df1 = %d, df2 = %.2f\n', df1_BF, df2_BF);
    
    F_critical_BF = finv(0.95, df1_BF, df2_BF);
    p_value_BF   = 1 - fcdf(F_BF, df1_BF, df2_BF);
    
    fprintf('F临界值(α=0.05) = %.6f\n', F_critical_BF);
    fprintf('p值 = %.6f\n', p_value_BF);
    
    fprintf('\n决策：\n');
    if F_BF > F_critical_BF
        fprintf('F = %.6f > %.6f，拒绝H0\n', F_BF, F_critical_BF);
        fprintf('结论：在α=0.05水平下，地区A和B的均值向量显著不同\n');
        BF_result = '显著不同';
    else
        fprintf('F = %.6f ≤ %.6f，接受H0\n', F_BF, F_critical_BF);
        fprintf('结论：在α=0.05水平下，地区A和B的均值向量无显著差异\n');
        BF_result = '无显著差异';
    end
else
    fprintf('\n样本量较小，f = %.2f ≤ p = %d\n', f_param, p);
    fprintf('使用渐近χ²分布方法：\n');
end