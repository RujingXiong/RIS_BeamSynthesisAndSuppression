% CreatTime: 2024.08.04
 % Complete:
 % Modified: 2024.09.18
 % E-mail: 
 % Description: Optimization on constrainted max-min problem.

clear all
clc
% 发射设置
TX.d = 5;
TX.phi = deg2rad(0);
TX.theta = deg2rad(0);
% 接收设置
% phi_out = deg2rad([0;180;180]);
% theta_out = deg2rad([20;20;50]); 
phi_out = deg2rad([0;180;180]);
theta_out = deg2rad([160;160;130]); 
d_out =[30;30;30];
% RIS设置
RIS_row = 1;
RIS_col = 16;
N = RIS_row*RIS_col;
frequency = 5.8e9;
% 约束设置
K = length(theta_out);
T = 1;
race = [1;2;4];
% 计算无约束情况下
Coeff.K = K;
Coeff.T = T;
Coeff.D = ones(K,1);
Coeff.N = N;
lb = zeros(N, 1);
ub = 2*pi*ones(N, 1);
Coeff.lb = lb;
Coeff.ub = ub;
% omega0 = lb + (ub-lb).*rand(N, 1);
omega0 = rand(1,1)*2*pi.*ones(N, 1);
h = cell(K,1);
for k=1:K
    Rx_position = d_out(k).*[sin(theta_out(k))*cos(phi_out(k));sin(theta_out(k))*sin(phi_out(k));cos(theta_out(k))];
    [h{k},R1,RIS_patch_position] = channel(TX, Rx_position, frequency, RIS_row, RIS_col); %Appropriate amplification is required during the computation process to facilitate MATLAB calculations
    h{k} = 1/race(k).*h{k}'*h{k};
end
Coeff.h = h;
[result, ~, Coeff] = TestMaxMin1(Coeff, omega0);
opt = exp(1i.*result.omega);
strength = cal(opt,TX,frequency,30, RIS_row, RIS_col);
strength = strength./(1e8); %recover the actual signal power
%strength = strength./max(strength); %normalize
xita = -90:0.25:90;
figure
plot(xita,strength);
grid on
xlabel('Angle({\circ})');
ylabel('Power');
grid on

%%
H = cell(T,1);
H{1} = Generate(R1, N, frequency, RIS_patch_position, 30, 0, 10, 179,181 );
sigma = 0.1*real(opt'*H{1}*opt); %This also requires appropriate amplification to facilitate the calculation.
%sigma = 0.001;
H{1} = (sigma/N).*eye(N,N) - H{1};

% H{2} = Generate(R1, N, frequency, RIS_patch_position, 30, 0, 15,179,181);
% sigma = 0.1*real(opt'*H{2}*opt);
% H{2} = (sigma/N).*eye(N,N) - H{2};
% 
% H{3} = Generate(R1, N, frequency, RIS_patch_position, 30, 122, 127);
% sigma = 0.2*real(opt'*H{3}*opt);
% H{3} = (sigma/N).*eye(N,N) - H{3};
Coeff.H = H;
[result, Coeff] = TestMaxMin(Coeff, omega0);
opt2 = exp(1i.*result.omega); % The optimal of constrained problem
strength = cal(opt2,TX,frequency,30, RIS_row, RIS_col);
strength = strength./(1e8);
%strength = strength./max(strength); %normalize
xita = -90:0.25:90;


figure
plot(xita,strength,'b-','LineWidth',1)
hold on
%xline(cos(TX.phi)*TX.theta,'r--','LineWidth',1)
xline(0,'r--','LineWidth',1);
xline(-10,'r--','LineWidth',1);
xlim([-90,90]);
%annotation('doublearrow', [0.45 0.5], [0.7 0.7]) % 绘制双向箭头([x1,x2],[y1,y2])
%legend('Incident','Reflect')
xlabel('Angle({\circ})')
ylabel('Power (mW)')
grid on

% polar dB figure
angle = (-90:0.25:90);
C1 = cal(opt2,TX,frequency,30, RIS_row, RIS_col);
C1 = C1./max(C1);
figure('NumberTitle','off','Name','Figure 11.9 in Polar','Position',[0 0 600 850]);
polardb(angle*pi/180,10*log10(abs(C1)),-30,'b');
hold on;
angle_of_interest = 0; % mark the incident angle
plot([angle_of_interest, angle_of_interest], [0, 3], 'r--','LineWidth',1.5); % use a red dashed line to highlight the incident angle, 表示从笛卡尔坐标系的(0,0)到(0,3)画一条线，
% 绘制从圆心指向 angle=10 的红虚线
plot([0, -0.5208], [0,2.954], 'r--', 'LineWidth', 1.5);%([y1,y2],[x1,x2])
hold off;
%title('Beampattern');
grid on;