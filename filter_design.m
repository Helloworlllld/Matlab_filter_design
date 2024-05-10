n_order = @(L_As, L_Ar, omega_s, omega_c) acosh(sqrt((10^(L_As/10) - 1) / (10^(L_Ar/10) - 1))) / acosh(omega_s / omega_c);

% 参数
L_As = 20; % 示例参数 L_As
L_Ar = 10; % 示例参数 L_Ar
omega_s = 2*pi*6e9; % 示例参数 omega_s
omega_c = 1; % 示例参数 omega_c

% 计算 n 的值
n_value = n_order(L_As, L_Ar, omega_s, omega_c);

% 显示结果
disp(['n >= ', num2str(n_value)]);