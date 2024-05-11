%切比雪夫Ⅰ型带通滤波器
%输入参数
w_center=2*pi*5e9;              %中心角频率
relative_bw = 0.1;              %相对带宽
w_stop = 2*pi*0.5e9;            %指定阻带偏离中心的角频率
Stop_dB = 15;                   %对应阻带角频率处的衰减
Reflect_dB = 20;                %通带内最小回波损耗
Z0 = 50;                        %特性阻抗


% 计算其余参数
w_pass = w_center*relative_bw;                                  %通带带宽角频率
Ripple_dB = abs(20*log10(sqrt(1-power(10,-2*Reflect_dB/20))));  %通带内插损纹波,取正数
epsilon = sqrt(power(10,(Ripple_dB/10)) - 1);                   %插损纹波对应的epsilon

% 阶数计算公式
n_order = @(L_As, L_Ar, w_stop, w_pass_half) acosh(sqrt((10^(L_As/10) - 1) / (10^(L_Ar/10) - 1))) / acosh(w_stop / w_pass_half);

% 计算阶数的值
n_value = n_order(Stop_dB, Ripple_dB, w_stop, w_pass/2);
% 向上取整
order = ceil(n_value);

% 显示阶数
disp(['order= ', num2str(order)]);

% 初始化归一化元件值g
g = zeros(1,order);

%阶数分奇偶
if mod(order,2) == 1
    % 奇数阶，直接通过公式计算g
    beta = log(coth(Ripple_dB/(40/log(10))));
    gamma = sinh(beta/(2*order));
    g(1) = (2/gamma)*sin(pi/(2*order));
    for i = 2:order
        g(i) = 1/g(i-1) * (4*sin((2*i-1)*pi/(2*order))*sin((2*i-3)*pi/(2*order))) / (gamma^2 + sin((i-1)*pi/order)^2);
    end
else
    % 偶数阶，为匹配阻抗，搬移切比雪夫多项式零点后使用定义直接计算阻抗，用辗转相除法求g
    syms s;%定义s为符号变量,s=1j*w(下面写作omega)
    %迭代法写出order阶切比雪夫多项式
    T_2 = 1;
    T_1 = s;
    T = 2*power(s,2) - 1;
    T_temp = 0;
    for i = 1:order-2
        T_temp = T;
        T = 2*s*T - T_1;
        T_2 = T_1;
        T_1 = T_temp;
    end
    T = collect(T,s); %整理多项式

    %搬移零点，保证直流增益为1
    omega = sqrt(-power(s,2));
    omega_a = cos((order-1)*pi/(2*order));
    omega_1 = sqrt(power(omega,2)*(1-power(omega_a,2))+power(omega_a,2));
    T = subs(T,s,omega_1);
    T_coeffs = coeffs(T,s);
    F_s = T*epsilon;    %F(s)定义式中包含P(s)，最后计算阻抗的比值式中会与E(s)中的P(s)约掉，故不考虑
    H_s = 1+power(epsilon,2)*power(T,2);%这里其实是H(s)的分母，分子为1，直接倒过来方便计算
    H_s = collect(H_s); %整理多项式
    
    % 寻找H(s)极点,即分母的零点
    original_den_sym = coeffs(H_s, s, "All");
    poles = roots(original_den_sym);
    poles = collect(poles);
    % 选择左半平面的极点
    left_half_poles = poles(real(poles) < 0);
    %digits(64); %解析转换数值，控制精度
    left_half_poles_vpa = vpa(left_half_poles);
    coeffs_vector = zeros(order,1);
    for i = 1:order
        coeffs_vector(i) = sym2poly(left_half_poles_vpa(i));
    end
    
    % 用左半平面极点构建传递函数的分母多项式，求系数
    new_den = poly(coeffs_vector);
    % 将系数根据最后一位归一化，保证输出负载阻抗为1
    new_den = 1/new_den(end)*new_den;
    %依据多项式系数写出新的传递函数多项式
    E_s = poly2sym(new_den, s);
    
    %统一F_s和E_s的精度(由于F(s)为syms，E(s)数值)
    F_s = vpa(F_s);
    E_s = vpa(E_s);
    
    %分情况写出阻抗分子分母多项式
    if mod(order,4) == 0
        A = E_s+F_s;
        B = E_s-F_s;
    else
        B = E_s+F_s;
        A = E_s-F_s;
    end
    %整理多项式
    A = collect(A);
    B = collect(B);
    
    %发现由于数值和解析值的转换误差，阶数大于2时B多出一项系数极小的order阶项，减去
    if order > 2
        coeffs_B = coeffs(B, s);
        B = B - coeffs_B(end)*power(s,order);
    end
    
    %辗转相除法求g
    g_s = sym('X',[1 order]);
    for i = 1:order
        [Q,R]= quorem(A, B, s);
        g_s(i) = Q;
        [A,B]=numden(B/R);
        gtemp = coeffs(g_s(i));
        g(i) = gtemp(end);
    end
end

% 计算L和C的值
L = zeros(1,order);
C = zeros(1,order);
for i = 1:order
    if mod(i,2) == 1
        C(i) = relative_bw/(w_center*Z0*g(i));
        L(i) = Z0*g(i)/(relative_bw*w_center);
    else
        L(i) = relative_bw*Z0/(w_center*g(i));
        C(i) = g(i)/(Z0*relative_bw*w_center);
    end
end

%逐组输出L和C的值
for i = 1:order
    disp(['L', num2str(i), ' = ', num2str(L(i)), ';']);
    disp(['C', num2str(i), ' = ', num2str(C(i)), ';']);
end
