%������ͼ�Ա�
clear;clc;
close all;
i = 1;%%�ڼ�������(1,2,0 �ֱ�����һ�εڶ��ε���������)
B=xlsread('Traj_20s_loop2_1',['Sheet',num2str(i)]);
In = 10;
q = B(1:In:end,1:6)';
dq = B(1:In:end,7:12)';
ddq = B(1:In:end,13:18)';
Current = B(1:In:end,37:42)';
CurrentFilt = B(1:In:end,43:48)';
Force = B(1:In:end,49:54)';
ForceFilt  = B(1:In:end,55:60)';
P1 = Force';
P2 = ForceFilt';
x_label = 'ʱ��t(0.001s)';
y_label = '����t(N)';
first_legend = '����ԭʼֵ';
second_legend = '�����˲�ֵ';
%% ��ͼ
%n = size(q,2);
n=1000;
A = figure(1);
set(A, 'unit', 'normalized', 'position', [0,0,1,1]);
subplot(2,3,1)
plot(1:n,P1(1:n,1),'b',1:n,P2(1:n,1),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('�ؽ�1');
hold on;
subplot(2,3,2)
plot(1:n,P1(1:n,2),'b',1:n,P2(1:n,2),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('�ؽ�2');
hold on;
subplot(2,3,3)
plot(1:n,P1(1:n,3),'b',1:n,P2(1:n,3),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('�ؽ�3');
hold on;
subplot(2,3,4)
plot(1:n,P1(1:n,4),'b',1:n,P2(1:n,4),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('�ؽ�4');
hold on;
subplot(2,3,5)
plot(1:n,P1(1:n,5),'b',1:n,P2(1:n,5),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('�ؽ�5');
hold on;
subplot(2,3,6)
plot(1:n,P1(1:n,6),'b',1:n,P2(1:n,6),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('�ؽ�6');
frame = getframe(A);