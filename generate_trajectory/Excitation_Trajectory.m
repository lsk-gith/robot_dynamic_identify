%激励轨迹设计与优化主函数
clear 
%% 自定义参数
joint_num = 6;  %关节数目
T = 20;    %轨迹运行周期T
dt = 0.01;   %采样间隔dt，dt=0.002导致数据点过多运行时间过久，因此设置dt=0.1
t = 0:dt:T; %时间t
N = size(t,2); %轨迹点数
w = 2*pi/T;    %频率w
Fourier_N = 3;  %傅里叶级数

%% 设置关节约束
% qmin:关节角度最小值; qmax:关节角度最大值; dqmin:关节角速度最大值; qmin:关节角加速度最大值;
qmin = [-120;-90;-120;-120;-120;-120]/180*pi;   %qmin = -ones(6,1)*120/180*pi;
qmax = [120;90;120;120;120;120]/180*pi; %qmax = ones(6,1)*120/180*pi;
dqmax = [100;80;80;80;80;80]/180*pi;%dqmax = ones(6,1)*80/180*pi; %dqmax = ones(6,1)*40/180*pi; 
ddqmax = [48;24;24;24;24;24]/180*pi;%ddqmax = ones(6,1)*24/180*pi; %ddqmax = ones(6,1)*12/180*pi; 
q_offset = (qmin+qmax)/2; %角度约束中点

%% 优化
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];

OPTIONS =optimoptions('fmincon','Algorithm','active-set');  

%设置傅里叶级数系数x的初始值
%1) 一开始,可以选择从随机数开始
%x_temp = -1+rand(joint_num*(2*Fourier_N+1),1)*2;
%2）从之前某次的优化结果开始
loadpath = './2020_1012/Traj_20s_3_temp1.mat';
x_temp =  load(loadpath);
x_temp = x_temp.x;


%优化目标为矩阵条件数COND(x)，约束为mycon(x), 优化傅里叶级数系数x
[x,fval,exitflag,output] = fmincon(@(x) COND(x),x_temp,[],[],[],[],lb,ub,@(x) mycon(x),OPTIONS);

%% 根据优化系数计算傅里叶级数激励轨迹
dt = 0.002;   %此时更改dt=0.002
t = 0:dt:T; %时间t
N = size(t,2); %轨迹点数

q = zeros(joint_num,N); dq = zeros(joint_num,N); ddq = zeros(joint_num,N);

for i = 1:joint_num
    for j = 1:Fourier_N
        q(i,:) = q(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))/(w*j)*sin(w*j*t) - x(2*j+(2*Fourier_N+1)*(i-1))/(w*j)*cos(w*j*t);
        dq(i,:) = dq(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))*cos(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*sin(w*j*t);
        ddq(i,:) = ddq(i,:) + (-x(2*j-1+(2*Fourier_N+1)*(i-1))*w*j*sin(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*w*j*cos(w*j*t));
    end
    q(i,:) = q(i,:) + x((2*Fourier_N+1)*i);
end

% 画图
figure,plot(t,q'*180/pi)
% figure,plot(dq'*180/pi)
% figure,plot(ddq'*180/pi)
%% 扩展成多周期
% q = [q,q]; dq = [dq,dq]; ddq = [ddq,ddq];
% t = [0:0.002:0.002*(size(q,2)-1)];