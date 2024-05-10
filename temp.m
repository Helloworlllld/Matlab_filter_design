% 定义原始函数的分母多项式
syms s;
H_w = 1+33.97*s^8+56.28*s^6+23.31*s^4;
original_den_sym = coeffs(H_w, s, "All");


% 计算极点
poles = roots(original_den_sym);

% 选择左半平面的极点
left_half_poles = poles(real(poles) < 0);
left_half_poles_double = double(left_half_poles);

% 进行部分分式展开（这里只是一个示例，具体的系数需要计算）
[r, p, k] = residue(1, left_half_poles_double);

% 构建新的传递函数的分母多项式
new_den = poly(left_half_poles_double);

