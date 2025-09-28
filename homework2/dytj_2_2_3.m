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
% 成分名称
components = {'SiO_2', 'Fe_2O_3', 'K_2O'};

%% 问题3：多个正态总体均值相等性检验（地区A、B、C）
fprintf('问题3：多个正态总体均值相等性检验\n');
fprintf('H0: μ1 = μ2 = μ3 vs H1: 至少有一对均值不相等\n\n');

% 基本参数
p = 3; 
r = 3;
n_total = n_A + n_B + n_C;

% 样本均值向量
X_bar_A = mean(data_A)';
X_bar_B = mean(data_B)';
X_bar_C = mean(data_C)';

% 总体均值向量
X_bar_total = (n_A * X_bar_A + n_B * X_bar_B + n_C * X_bar_C) / n_total;

fprintf('各组样本统计量：\n');
fprintf('地区A: n1 = %d, 均值向量 = [%.4f, %.4f, %.4f]\n', n_A, X_bar_A(1), X_bar_A(2), X_bar_A(3));
fprintf('地区B: n2 = %d, 均值向量 = [%.4f, %.4f, %.4f]\n', n_B, X_bar_B(1), X_bar_B(2), X_bar_B(3));
fprintf('地区C: n3 = %d, 均值向量 = [%.4f, %.4f, %.4f]\n', n_C, X_bar_C(1), X_bar_C(2), X_bar_C(3));
fprintf('总体均值估计: X̄ = [%.4f, %.4f, %.4f]\n\n', X_bar_total(1), X_bar_total(2), X_bar_total(3));

% 将数据重新组织成标准格式
% X^(k)_j 表示第k个总体的第j个观测向量
groups = {data_A, data_B, data_C};
n_groups = [n_A, n_B, n_C];
X_bars = {X_bar_A, X_bar_B, X_bar_C};

% 计算B = SS(TR) - 组间平方和矩阵
% B = Σ(k=1 to r) n_k(X̄_k - X̄)(X̄_k - X̄)'
fprintf('\n组间平方和矩阵 B = SS(TR)：\n');
B = zeros(p, p);
for k = 1:length(groups)
    diff_k = X_bars{k} - X_bar_total;
    B_k = n_groups(k) * (diff_k * diff_k');
    B = B + B_k;
end
fprintf('B = SS(TR) = \n');
disp(B);

% 计算E = SSE - 组内平方和矩阵（误差平方和）
% E = Σ(k=1 to r)Σ(j=1 to n_k)(X^(k)_j - X̄_k)(X^(k)_j - X̄_k)'
fprintf('\n组内平方和矩阵 E = SSE：\n');
E = zeros(p, p);
for k = 1:length(groups)
    E_k = zeros(p, p);
    data_k = groups{k};
    X_bar_k = X_bars{k};
    
    for j = 1:n_groups(k)
        diff_kj = data_k(j, :)' - X_bar_k;
        E_k = E_k + diff_kj * diff_kj';
    end
    E = E + E_k;
end
fprintf('E = SSE = \n');
disp(E);

% 计算W = SST - 总平方和矩阵
% W = Σ(k=1 to r)Σ(j=1 to n_k)(X^(k)_j - X̄)(X^(k)_j - X̄)'
fprintf('\n总平方和矩阵 W = SST：\n');
W = zeros(p, p);
for k = 1:length(groups)
    data_k = groups{k};
    for j = 1:n_groups(k)
        diff_total = data_k(j, :)' - X_bar_total;
        W = W + diff_total * diff_total';
    end
end
fprintf('W = SST = \n');
disp(W);

%% Wilks Lambda(Λ) 检验
fprintf('Wilks Lambda (Λ) 检验：\n');
Lambda = det(E) / det(W);
fprintf('   Λ = |E|/|W| = %.8f\n', Lambda);

% 使用Bartlett近似转换为F分布
% 当p和r-1都较小时的精确F分布转换
if p == 1 || r == 2
    % 精确F分布
    F_Lambda = ((n_total - r - p + 1)/p) * ((1 - Lambda)/Lambda);
    df1_Lambda = p;
    df2_Lambda = N - r - p + 1;
    fprintf('   使用精确F分布\n');
else
    % 使用Rao's R统计量近似
    w = n_total - r - (p - r + 2)/2;
    t = sqrt((p^2*(r-1)^2 - 4)/(p^2 + (r-1)^2 - 5));
    if p^2 + (r-1)^2 - 5 <= 0, t = 1; end
    
    df1_Lambda = p * (r - 1);
    df2_Lambda = w * t - (p*(r-1) - 2)/2;
    
    F_Lambda = ((1 - Lambda^(1/t))/(Lambda^(1/t))) * (df2_Lambda/df1_Lambda);
    fprintf('   使用Raos R近似，t = %.4f\n', t);
end

fprintf('   F统计量 = %.6f\n', F_Lambda);
fprintf('   自由度: df1 = %d, df2 = %.2f\n', df1_Lambda, df2_Lambda);

F_crit_Lambda = finv(0.95, df1_Lambda, df2_Lambda);
p_val_Lambda = 1 - fcdf(F_Lambda, df1_Lambda, df2_Lambda);

fprintf('   F临界值(α=0.05) = %.6f\n', F_crit_Lambda);
fprintf('   p值 = %.6f\n', p_val_Lambda);

if F_Lambda > F_crit_Lambda
    fprintf('   \n决策: \n   F = %.4f > %.4f，拒绝H0\n', F_Lambda, F_crit_Lambda);
    fprintf('   结论: 三个地区的均值向量不全相等\n\n');
    result_Lambda = '拒绝H0';
else
    fprintf('   \n决策: \n   F = %.4f ≤ %.4f，接受H0\n', F_Lambda, F_crit_Lambda);
    fprintf('   结论: 三个地区的均值向量相等\n\n');
    result_Lambda = '接受H0';
end