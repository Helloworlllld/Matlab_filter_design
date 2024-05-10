syms s; %参数s=jw

% 定义电感和电容的值
L1 = 39.73e-9;
C1 = 159.4e-15;
C2 = 289.5e-12;
L2 = 21.88e-12;
C3 = 87.05e-15;
L3 = 72.75e-9;
C4 = 289.5e-12;
L4 = 21.88e-12;
C5 = 159.4e-15;
L5 = 39.73e-9;

%系统为50Ω阻抗匹配
%总输入阻抗，计入激励源50Ω内阻和输出50Ω负载
Zin = 50+s*L1+(1/(s*C1))+1/((1/(s*C2))+s*L2+1/((1/(s*C3))+s*L3+1/((1/(s*C4))+s*L4+(1/(50+s*L5+(1/(s*C5)))))));

%用阻抗模型计算输出电压,引入等效电流,同时将Ui计为1(若不为1最后也会约掉)方便计算
I_total=2/Zin;
ZP1=1/(1/(s*C2)+s*L2);
ZP2=1/(s*C3)+s*L3+1/((1/(s*C4))+s*L4+(1/(50+s*L5+(1/(s*C5)))));
I_1=I_total/(1+ZP2/ZP1);
ZP3=1/(1/(s*C4)+s*L4);
ZP4=50+s*L5+(1/(s*C5));
I_2=I_1/(1+ZP4/ZP3);
H=I_2*50;


fmin = 1;     % 最小频率1Hz
fmax = 5e9;   % 最大频率5GHz
npoints = 1000; % 频率点数
f = linspace(fmin,fmax,npoints); % 线性频率向量
w = 2*pi*f;   % 对应的角频率向量
H_jw = zeros(size(w));
reflect = zeros(size(w));
S1221_dB = zeros(size(w));
S11_dB = zeros(size(w));
for i = 1:length(w)
   H_jw(i) = abs(subs(H,s,1j*w(i))); 
   reflect(i) = abs(sqrt(1-power(H_jw(i), 2)));
   S1221_dB(i) = 20*log10(H_jw(i));
   S11_dB(i) = 20*log10(reflect(i));
end


figure;
plot(f/1e9, S11_dB, 'LineWidth', 2, 'DisplayName', 'S11');
hold on;
plot(f/1e9, S1221_dB, 'LineWidth', 2, 'DisplayName', 'S12/S21');
grid on;
xlim([fmin/1e9 fmax/1e9]);
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('S11 and S12/21 (0-5GHz)');
legend;