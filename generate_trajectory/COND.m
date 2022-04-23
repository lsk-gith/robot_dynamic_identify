function [cond_val] = COND(x)
%COND2(x) 返回轨迹的所有信息矩阵的集合K_total的条件数
%x为傅里叶级数系数，大小为 joint_num*(2*Fourier_N+1)x1
%x = [a11; b11; a12; b12; ... q10; ...]
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

%% 计算回归矩阵条件数
% modified_dh参数
alpha=[0;pi/2;0;0;pi/2;-pi/2];    
a=[0;0;-0.264;-0.237;0;0];         
d=[0.144;0;0;0.1065;0.114;0.09];           
theta=[0;-pi/2;0;-pi/2;0;0];
dh_list = [alpha a d theta];
%重力加速度g
g = 9.81;
%摩擦力项数
nf = 2;
%计算激轨迹
q = zeros(joint_num,N); dq = zeros(joint_num,N); ddq = zeros(joint_num,N);
for i = 1:joint_num
    for j = 1:Fourier_N
        q(i,:) = q(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))/(w*j)*sin(w*j*t) - x(2*j+(2*Fourier_N+1)*(i-1))/(w*j)*cos(w*j*t);
        dq(i,:) = dq(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))*cos(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*sin(w*j*t);
        ddq(i,:) = ddq(i,:) + (-x(2*j-1+(2*Fourier_N+1)*(i-1))*w*j*sin(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*w*j*cos(w*j*t));
    end
    q(i,:) = q(i,:) + x((2*Fourier_N+1)*i);
end
%根据激励轨迹计算回归矩阵
J = zeros(joint_num,6,N); K = zeros(joint_num,10*joint_num,N); Kf = zeros(joint_num,2*joint_num,N);
for cnt = 1:N
[J(:,:,cnt),K(:,:,cnt),Kf(:,:,cnt)] = Compute_Dynmatrix(q(:,cnt), dq(:,cnt), ddq(:,cnt), g, dh_list,nf);
end
%K_total为所有时刻回归矩阵的集合
K_total = zeros(joint_num*N,10*joint_num);
for i = 1:N
    K_total(joint_num*(i-1)+1:joint_num*i,:) = K(:,:,i); 
end

%删去K_total中一直为0的列，1~9，11，14，否则条件数为无穷inf
K_total = [K_total(:,10), K_total(:,12:13), K_total(:,15:end)];
%K_total1 = K_total(:,find(sum(abs(K_total))'>sqrt(eps)));
%进一步qr分解，得最小参数集回归矩阵KK
[~ , ~, P] = qr(K_total);
nd = 36; %nd = rank(K_total); 
KK_total = K_total*P;
KK = KK_total(:,1:nd);
%计算最小参数集回归矩阵KK的条件数
cond_val =cond(KK)

%% 优化过程记录合适的傅里叶级数系数x
% 若有满足关节约束且条件数小于一定值(如小于100), 则将傅里叶系数等值保存至给定路径savepath下的mat数组
% 若满足关节约束, 则将最小回归矩阵条件数值和对应傅里叶系数写入给定路径savepath2的txt文件中
% 优化超时还找不到满足条件数小于指定值的傅里叶系数，则需要手动在txt文件中挑选合适的傅里叶系数, 计算傅里叶级数激励轨迹。

cond_max_val = 100;
savepath = './2020_1012/Traj_fuhe_001_20s_3_newcon1_80.mat';
savepath2 = './2020_1012/cond_fuhe_001_20s_3_newnewcon.txt';

if((abs(q(:,1)*180/pi)<1).*(abs(dq(:,1)*180/pi)<1).*(abs(ddq(:,1)*180/pi)<1)...  %初始关节角度小于1度(即初始关节角度为0)
    .*(qmin<min(q,[],2)).*(max(q,[],2) < qmax).*(max(abs(dq),[],2) < dqmax).*(max(abs(ddq),[],2) < ddqmax))
    fid = fopen(savepath2,'a+');
    fprintf(fid,'%f',cond_val);
    for i = 1 : length(x)
        fprintf(fid,',%f',x(i));
    end
    fprintf(fid,'\n');
    fclose(fid);

    if(cond_val < cond_max_val)
        save(savepath,'q','dq','ddq','x','t');
        ni=input('条件数小于100，继续请输入一个数字，并按回车：');
        pause;
    end
end
end