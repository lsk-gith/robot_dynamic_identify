%�����켣������Ż�������
clear 
%% �Զ������
joint_num = 6;  %�ؽ���Ŀ
T = 20;    %�켣��������T
dt = 0.01;   %�������dt��dt=0.002�������ݵ��������ʱ����ã��������dt=0.1
t = 0:dt:T; %ʱ��t
N = size(t,2); %�켣����
w = 2*pi/T;    %Ƶ��w
Fourier_N = 3;  %����Ҷ����

%% ���ùؽ�Լ��
% qmin:�ؽڽǶ���Сֵ; qmax:�ؽڽǶ����ֵ; dqmin:�ؽڽ��ٶ����ֵ; qmin:�ؽڽǼ��ٶ����ֵ;
qmin = [-120;-90;-120;-120;-120;-120]/180*pi;   %qmin = -ones(6,1)*120/180*pi;
qmax = [120;90;120;120;120;120]/180*pi; %qmax = ones(6,1)*120/180*pi;
dqmax = [100;80;80;80;80;80]/180*pi;%dqmax = ones(6,1)*80/180*pi; %dqmax = ones(6,1)*40/180*pi; 
ddqmax = [48;24;24;24;24;24]/180*pi;%ddqmax = ones(6,1)*24/180*pi; %ddqmax = ones(6,1)*12/180*pi; 
q_offset = (qmin+qmax)/2; %�Ƕ�Լ���е�

%% �Ż�
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];

OPTIONS =optimoptions('fmincon','Algorithm','active-set');  

%���ø���Ҷ����ϵ��x�ĳ�ʼֵ
%1) һ��ʼ,����ѡ����������ʼ
%x_temp = -1+rand(joint_num*(2*Fourier_N+1),1)*2;
%2����֮ǰĳ�ε��Ż������ʼ
loadpath = './2020_1012/Traj_20s_3_temp1.mat';
x_temp =  load(loadpath);
x_temp = x_temp.x;


%�Ż�Ŀ��Ϊ����������COND(x)��Լ��Ϊmycon(x), �Ż�����Ҷ����ϵ��x
[x,fval,exitflag,output] = fmincon(@(x) COND(x),x_temp,[],[],[],[],lb,ub,@(x) mycon(x),OPTIONS);

%% �����Ż�ϵ�����㸵��Ҷ���������켣
dt = 0.002;   %��ʱ����dt=0.002
t = 0:dt:T; %ʱ��t
N = size(t,2); %�켣����

q = zeros(joint_num,N); dq = zeros(joint_num,N); ddq = zeros(joint_num,N);

for i = 1:joint_num
    for j = 1:Fourier_N
        q(i,:) = q(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))/(w*j)*sin(w*j*t) - x(2*j+(2*Fourier_N+1)*(i-1))/(w*j)*cos(w*j*t);
        dq(i,:) = dq(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))*cos(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*sin(w*j*t);
        ddq(i,:) = ddq(i,:) + (-x(2*j-1+(2*Fourier_N+1)*(i-1))*w*j*sin(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*w*j*cos(w*j*t));
    end
    q(i,:) = q(i,:) + x((2*Fourier_N+1)*i);
end

% ��ͼ
figure,plot(t,q'*180/pi)
% figure,plot(dq'*180/pi)
% figure,plot(ddq'*180/pi)
%% ��չ�ɶ�����
% q = [q,q]; dq = [dq,dq]; ddq = [ddq,ddq];
% t = [0:0.002:0.002*(size(q,2)-1)];