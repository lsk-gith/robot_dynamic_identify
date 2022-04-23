function [c,ceq] = mycon(x)
%设置不等式约束c和等式约束ceq, 由推导无需求取具体轨迹
%x为傅里叶级数系数，大小为 6x(2*Fourier_N+1)
%x = [a11 b11 a12 b12 ... q10; ...]
%% 自定义参数(同 Excitation_Trajectory.m)
joint_num = 6;  %关节数目
T = 20;    %轨迹运行周期T
dt = 0.01;   %采样间隔dt，dt=0.002导致数据点过多运行时间过久，因此设置dt=0.1
t = 0:dt:T; %时间t
N = size(t,2); %轨迹点数
w = 2*pi/T;    %频率w
Fourier_N = 3;  %傅里叶级数

%% 设置关节约束(同 Excitation_Trajectory.m)
% qmin:关节角度最小值; qmax:关节角度最大值; dqmin:关节角速度最大值; qmin:关节角加速度最大值;
qmin = [-120;-90;-120;-120;-120;-120]/180*pi;   %qmin = -ones(6,1)*120/180*pi;
qmax = [120;90;120;120;120;120]/180*pi; %qmax = ones(6,1)*120/180*pi;
dqmax = [100;80;80;80;80;80]/180*pi;%dqmax = ones(6,1)*80/180*pi; %dqmax = ones(6,1)*40/180*pi; 
ddqmax = [48;24;24;24;24;24]/180*pi;%ddqmax = ones(6,1)*24/180*pi; %ddqmax = ones(6,1)*12/180*pi; 
q_offset = (qmin+qmax)/2; %角度约束中点

%% 计算激励轨迹轨迹
q = zeros(joint_num,N); %dq = zeros(joint_num,N); ddq = zeros(joint_num,N);

for i = 1:joint_num
    for j = 1:Fourier_N
        q(i,:) = q(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))/(w*j)*sin(w*j*t) - x(2*j+(2*Fourier_N+1)*(i-1))/(w*j)*cos(w*j*t);
%        dq(i,:) = dq(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))*cos(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*sin(w*j*t);
%        ddq(i,:) = ddq(i,:) + (-x(2*j-1+(2*Fourier_N+1)*(i-1))*w*j*sin(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*w*j*cos(w*j*t));
    end
    q(i,:) = q(i,:) + x((2*Fourier_N+1)*i);
end
%% 设置不等式约束和等式约束，保证角度、角速度和角加速度约束
c1 = zeros(joint_num,1); c2 = zeros(joint_num,1); c3 = zeros(joint_num,1);
ceq1 = zeros(joint_num,1); ceq2 = zeros(joint_num,1); ceq3 = zeros(joint_num,1);

for i = 1:joint_num
    for j = 1:Fourier_N
       c1(i) = c1(i)+1/(j*w)*sqrt(x(2*j-1+(2*Fourier_N+1)*(i-1)).^2+x(2*j+(2*Fourier_N+1)*(i-1)).^2) ;
       c2(i) = c2(i)+sqrt(x(2*j-1+(2*Fourier_N+1)*(i-1)).^2+x(2*j+(2*Fourier_N+1)*(i-1)).^2);
       c3(i) = c3(i)+w*j*sqrt(x(2*j-1+(2*Fourier_N+1)*(i-1)).^2+x(2*j+(2*Fourier_N+1)*(i-1)).^2);
       
       ceq1(i) = ceq1(i)+x(2*j+(2*Fourier_N+1)*(i-1))/(w*j);
       ceq2(i) = ceq2(i)+x(2*j-1+(2*Fourier_N+1)*(i-1));
       ceq3(i) = ceq3(i)+w*j*x(2*j+(2*Fourier_N+1)*(i-1));
    end
    c1(i) = c1(i)+abs(x((2*Fourier_N+1)*i)-q_offset(i))- (qmax(i)-q_offset(i));
    c2(i) = c2(i)- dqmax(i);
    c3(i) = c3(i)- ddqmax(i);
    
    ceq1(i) = ceq1(i)-x((2*Fourier_N+1)*i);
end

c = [c1;c2;c3];
ceq = [ceq1;ceq2;ceq3];
end