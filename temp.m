syms s;

% 定义原始函数的分母多项式
H_w = 1/sqrt(1+33.97*s^8+56.28*s^6+23.31*s^4);
original_den_sym = coeffs(1/H_w^2, s, "All");

% 计算极点
poles = roots(original_den_sym);

% 选择左半平面的极点
left_half_poles = poles(real(poles) < 0);
left_half_poles_vpa = vpa(left_half_poles);

% 使用 sym2poly 函数将符号表达式转换为系数向量
coeffs_vector = zeros(4,1);
for i = 1:4
    coeffs_vector(i) = sym2poly(left_half_poles_vpa(i));
end
new_den = poly(coeffs_vector)
% 构建新的传递函数的分母多项式
new_den = 1/new_den(end)*new_den;
