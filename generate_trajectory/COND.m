function [cond_val] = COND(x)
%COND2(x) ���ع켣��������Ϣ����ļ���K_total��������
%xΪ����Ҷ����ϵ������СΪ joint_num*(2*Fourier_N+1)x1
%x = [a11; b11; a12; b12; ... q10; ...]
%% �Զ������(ͬ Excitation_Trajectory.m)
joint_num = 6;  %�ؽ���Ŀ
T = 20;    %�켣��������T
dt = 0.01;   %�������dt��dt=0.002�������ݵ��������ʱ����ã��������dt=0.1
t = 0:dt:T; %ʱ��t
N = size(t,2); %�켣����
w = 2*pi/T;    %Ƶ��w
Fourier_N = 3;  %����Ҷ����

%% ���ùؽ�Լ��(ͬ Excitation_Trajectory.m)
% qmin:�ؽڽǶ���Сֵ; qmax:�ؽڽǶ����ֵ; dqmin:�ؽڽ��ٶ����ֵ; qmin:�ؽڽǼ��ٶ����ֵ;
qmin = [-120;-90;-120;-120;-120;-120]/180*pi;   %qmin = -ones(6,1)*120/180*pi;
qmax = [120;90;120;120;120;120]/180*pi; %qmax = ones(6,1)*120/180*pi;
dqmax = [100;80;80;80;80;80]/180*pi;%dqmax = ones(6,1)*80/180*pi; %dqmax = ones(6,1)*40/180*pi; 
ddqmax = [48;24;24;24;24;24]/180*pi;%ddqmax = ones(6,1)*24/180*pi; %ddqmax = ones(6,1)*12/180*pi; 
q_offset = (qmin+qmax)/2; %�Ƕ�Լ���е�

%% ����ع����������
% modified_dh����
alpha=[0;pi/2;0;0;pi/2;-pi/2];    
a=[0;0;-0.264;-0.237;0;0];         
d=[0.144;0;0;0.1065;0.114;0.09];           
theta=[0;-pi/2;0;-pi/2;0;0];
dh_list = [alpha a d theta];
%�������ٶ�g
g = 9.81;
%Ħ��������
nf = 2;
%���㼤�켣
q = zeros(joint_num,N); dq = zeros(joint_num,N); ddq = zeros(joint_num,N);
for i = 1:joint_num
    for j = 1:Fourier_N
        q(i,:) = q(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))/(w*j)*sin(w*j*t) - x(2*j+(2*Fourier_N+1)*(i-1))/(w*j)*cos(w*j*t);
        dq(i,:) = dq(i,:) + x(2*j-1+(2*Fourier_N+1)*(i-1))*cos(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*sin(w*j*t);
        ddq(i,:) = ddq(i,:) + (-x(2*j-1+(2*Fourier_N+1)*(i-1))*w*j*sin(w*j*t) + x(2*j+(2*Fourier_N+1)*(i-1))*w*j*cos(w*j*t));
    end
    q(i,:) = q(i,:) + x((2*Fourier_N+1)*i);
end
%���ݼ����켣����ع����
J = zeros(joint_num,6,N); K = zeros(joint_num,10*joint_num,N); Kf = zeros(joint_num,2*joint_num,N);
for cnt = 1:N
[J(:,:,cnt),K(:,:,cnt),Kf(:,:,cnt)] = Compute_Dynmatrix(q(:,cnt), dq(:,cnt), ddq(:,cnt), g, dh_list,nf);
end
%K_totalΪ����ʱ�̻ع����ļ���
K_total = zeros(joint_num*N,10*joint_num);
for i = 1:N
    K_total(joint_num*(i-1)+1:joint_num*i,:) = K(:,:,i); 
end

%ɾȥK_total��һֱΪ0���У�1~9��11��14������������Ϊ����inf
K_total = [K_total(:,10), K_total(:,12:13), K_total(:,15:end)];
%K_total1 = K_total(:,find(sum(abs(K_total))'>sqrt(eps)));
%��һ��qr�ֽ⣬����С�������ع����KK
[~ , ~, P] = qr(K_total);
nd = 36; %nd = rank(K_total); 
KK_total = K_total*P;
KK = KK_total(:,1:nd);
%������С�������ع����KK��������
cond_val =cond(KK)

%% �Ż����̼�¼���ʵĸ���Ҷ����ϵ��x
% ��������ؽ�Լ����������С��һ��ֵ(��С��100), �򽫸���Ҷϵ����ֵ����������·��savepath�µ�mat����
% ������ؽ�Լ��, ����С�ع����������ֵ�Ͷ�Ӧ����Ҷϵ��д�����·��savepath2��txt�ļ���
% �Ż���ʱ���Ҳ�������������С��ָ��ֵ�ĸ���Ҷϵ��������Ҫ�ֶ���txt�ļ�����ѡ���ʵĸ���Ҷϵ��, ���㸵��Ҷ���������켣��

cond_max_val = 100;
savepath = './2020_1012/Traj_fuhe_001_20s_3_newcon1_80.mat';
savepath2 = './2020_1012/cond_fuhe_001_20s_3_newnewcon.txt';

if((abs(q(:,1)*180/pi)<1).*(abs(dq(:,1)*180/pi)<1).*(abs(ddq(:,1)*180/pi)<1)...  %��ʼ�ؽڽǶ�С��1��(����ʼ�ؽڽǶ�Ϊ0)
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
        ni=input('������С��100������������һ�����֣������س���');
        pause;
    end
end
end