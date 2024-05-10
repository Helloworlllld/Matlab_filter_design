% 参数
w_center=2*pi*6e9;              %中心角频率
relative_bw = 0.1;              %相对带宽
w_pass = w_center*relative_bw;  %通带带宽角频率
w_stop = 2*pi*0.5e9;            %指定阻带偏离中心的角频率
Stop_dB = 15;                   %对应阻带角频率处的衰减
Reflect_dB = 20;                %通带内最小回波损耗
Ripple_dB = abs(20*log10(sqrt(1-power(10,-2*Reflect_dB/20))));%通带内插损纹波,取正数

%阶数计算公式
n_order = @(L_As, L_Ar, w_stop, w_pass_half) acosh(sqrt((10^(L_As/10) - 1) / (10^(L_Ar/10) - 1))) / acosh(w_stop / w_pass_half);

% 计算阶数的值
n_value = n_order(Stop_dB, Ripple_dB, w_stop, w_pass/2);
order = ceil(n_value);

% 显示阶数
disp(['order= ', num2str(order)]);

g = zeros(1,order+1);

% 计算g的值
%\beta=\ln\biggl[\coth\biggl(\frac{L_{\text{Ar}}}{17.37}\biggr)\biggr]
beta = log(coth(Ripple_dB/17.37));

%\gamma=\sinh\left(\frac\beta{2n}\right)
gamma = sinh(beta/(2*order));

g(1) = (2/gamma)*sin(pi/(2*order));

for i = 2:order
    g(i) = 1/g(i-1) * (4*sin((2*i-1)*pi/(2*order))*sin((2*i-3)*pi/(2*order))) / (gamma^2 + sin((i-1)*pi/order)^2);
end

%g_{n+1}=\begin{cases}1.0&\text{,}n\text{为奇数}\\\coth^2\biggl(\frac{\beta}{4}\biggr)&\text{,}n\text{为偶数}\end{cases}
if mod(order,2) == 1
    g(order+1) = 1.0;
else
    g(order+1) = coth(beta/4)^2;
end

