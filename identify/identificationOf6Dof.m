%�˺����ǲ�����ʶ
%% ����׼��
clear;clc;%%����ռ�
alpha=[0,pi/2,0,0,pi/2,-pi/2];    
a=[0,0,-0.264,-0.237,0,0];  
d=[0.144,0,0,0.1065,0.114,0.0895]; 
theta=[0,-90,0,-90,0,0]/180*pi;
dh=[alpha;a;d;theta]';
n = size(dh,1);

%% �������ݴ���
i = 1;%%�ڼ�������(1,2,0 �ֱ�����һ�εڶ��ε���������)
A=xlsread('Traj_20s_loop2_1',['Sheet',num2str(i)]);
%%%���ֵΪq,dq,ddq,forcefilt,force,currentfilt,current,
q = A(1:10:end,1:6)';%�������е�����Ӧ����N*6����ʽ
dq = A(1:10:end,7:12)';
ddq = A(1:10:end,13:18)';
CurrentInitial = A(1:10:end,37:42)';
Current = A(1:10:end,43:48)';
Force = A(1:10:end,55:60)';
tau = [];
for j=1:size(Current,2)
    tau = [tau;Current(:,j)];
end
%% ����KK
[KK,F] = calculteBigK(dh,q,dq,ddq);
%% ��ʶBP�ǽ��,Bp��sym��ʽ����С���Բ�����
[op1,op2,BP,byta,P]=qrLeast(KK,tau,F);
%% ������С���Բ�������ĸ��ʽ�Ľ�
[Bp,ParameterStand] = calculteBp(op1,op2,byta);
[G,M]=calculteGandM(dh,BP,op2,P);
[CC]=calculteCC(dh,G,op2,P,BP);
[BB]=calculteBB(dh,G,P,CC,BP,op2);
[C] = mergeBandC(CC,BB);
%% ��������������֤
i = 1;%%�ڼ�������(1,2,0 �ֱ�����һ�εڶ��ε���������)
B=xlsread('Traj_15s_loop3_2',['Sheet',num2str(i)]);
In = 10;
q = B(1:In:end,1:6)';
dq = B(1:In:end,7:12)';
ddq = B(1:In:end,13:18)';
Current = B(1:In:end,37:42)';
CurrentFilt = B(1:In:end,43:48)';
Force = B(1:In:end,49:54)';
ForceFilt  = B(1:In:end,55:60)';
tau = [];
for j=1:size(Current,2)
    tau = [tau;Current(:,j)];
end
%% 
[KK,FF] = calculteBigK(dh,q,dq,ddq);
%%% ��ʶ
kk = [];
for i = 1:size(op2,2)
    kk = [kk,KK(:,op2(i))];
end
% [Q,R,PP]=qr(kk);
kk = kk * P;
Ff = FF * BP(size(op2,2)+1:end);
tauIden = kk * BP(1:size(op2,2)) + Ff;
Fram = plotverify(tau,tauIden);
S = calculateS(tau,tauIden);
im1=frame2im(Fram);
% imwrite(im1,'identification4.png','png');

%% ʹ�÷���M G C ����֤
% tauIdentify = [];
% for i = 1:size(q,2)
% %     tauIdentify = [tauIdentify;double(subs(M,[sym('q',[6,1]);sym('dq',[6,1]);sym('Bp',[size(op2,2),1])],[q(:,i);dq(:,i);BP(1:size(op2,2))]))*ddq(:,i)+double(subs(G,[sym('q',[6,1]);sym('Bp',[size(op2,2),1])],[q(:,i);BP(1:size(op2,2))]))+ double(subs(C,[sym('q',[6,1]);sym('dq',[6,1]);sym('Bp',[size(op2,2),1])],[q(:,i);dq(:,i);BP(1:size(op2,2))]))*dq(:,i)];
%     tauIdentify = [tauIdentify;double(subs(M,[sym('q',[6,1]);sym('dq',[6,1]);sym('Bp',[size(op2,2),1])],[q(:,i);dq(:,i);BP(1:size(op2,2))]))*ddq(:,i)+ ...
%         double(subs(G,[sym('q',[6,1]);sym('Bp',[size(op2,2),1])],[q(:,i);BP(1:size(op2,2))]))+ ...
%         double(subs(CC,[sym('q',[6,1]);sym('Bp',[size(op2,2),1])],[q(:,i);BP(1:size(op2,2))]))*[dq(1,i)*dq(1,i);dq(2,i)*dq(2,i);dq(3,i)*dq(3,i);dq(4,i)*dq(4,i);dq(5,i)*dq(5,i);dq(6,i)*dq(6,i)]+ ...
%         double(subs(BB,[sym('q',[6,1]);sym('Bp',[size(op2,2),1])],[q(:,i);BP(1:size(op2,2))]))*[dq(1,i)*dq(2,i);dq(1,i)*dq(3,i);dq(1,i)*dq(4,i);dq(1,i)*dq(5,i);dq(1,i)*dq(6,i); ...
%         dq(2,i)*dq(3,i);dq(2,i)*dq(4,i);dq(2,i)*dq(5,i);dq(2,i)*dq(6,i); ...
%         dq(3,i)*dq(4,i);dq(3,i)*dq(5,i);dq(3,i)*dq(6,i); ...
%         dq(4,i)*dq(5,i); dq(4,i)*dq(6,i);...
%         dq(5,i)*dq(6,i)]];
%     i
% end
% tauIdentify = tauIdentify + Ff;
% Fram = plotverify(tau,tauIdentify);
% im2=frame2im(Fram);
% imwrite(im2,'identification5.png','png');


