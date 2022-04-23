function [frame] = plotverify(taulist1,taulist2)
m = size(taulist1,1);
n = m/6;
for i = 1:6
    for j = 1:m/6
        tau_test(j,i) = taulist1(6*(j-1)+i); 
        tau_real(j,i) = taulist2(6*(j-1)+i);
    end
end



P1 = tau_test;
P2 = tau_real;
x_label = '时间t(0.001s)';
y_label = '电流I(A)';
first_legend = '电流原始值';
second_legend = '电流辨识值';
%% 画图
%n = size(q,2);
n=800;
A = figure(1);
set(A, 'unit', 'normalized', 'position', [0,0,1,1]);
subplot(2,3,1)
plot(1:n,P1(1:n,1),'b',1:n,P2(1:n,1),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('关节1');
hold on;
subplot(2,3,2)
plot(1:n,P1(1:n,2),'b',1:n,P2(1:n,2),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('关节2');
hold on;
subplot(2,3,3)
plot(1:n,P1(1:n,3),'b',1:n,P2(1:n,3),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('关节3');
hold on;
subplot(2,3,4)
plot(1:n,P1(1:n,4),'b',1:n,P2(1:n,4),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('关节4');
hold on;
subplot(2,3,5)
plot(1:n,P1(1:n,5),'b',1:n,P2(1:n,5),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('关节5');
hold on;
subplot(2,3,6)
plot(1:n,P1(1:n,6),'b',1:n,P2(1:n,6),'r');
xlabel(x_label);
ylabel(y_label);
legend(first_legend,second_legend);
title('关节6');
frame = getframe(A);



% %% 力矩对比
% A = figure(1);
% set(A, 'unit', 'normalized', 'position', [0,0,1,1]);
% subplot(2,3,1)
% plot(1:n,tau_test(:,1),'b',1:n,tau_real(:,1),'r');
% xlabel('time（s）');
% ylabel('torque（N/m）');
% legend('tautest1','taulist1')
% hold on;
% subplot(2,3,2)
% plot(1:n,tau_test(:,2),'b',1:n,tau_real(:,2),'r');
% xlabel('time（s）');
% ylabel('torque（N/m）');
% legend('tautest2','taulist2')
% hold on;
% subplot(2,3,3)
% plot(1:n,tau_test(:,3),'b',1:n,tau_real(:,3),'r');
% xlabel('time（s）');
% ylabel('torque（N/m）');
% legend('tautest3','taulist3')
% hold on;
% subplot(2,3,4)
% plot(1:n,tau_test(:,4),'b',1:n,tau_real(:,4),'r');
% xlabel('time（s）');
% ylabel('torque（N/m）');
% legend('tautest4','taulist4')
% hold on;
% subplot(2,3,5)
% plot(1:n,tau_test(:,5),'b',1:n,tau_real(:,5),'r');
% xlabel('time（s）');
% ylabel('torque（N/m）');
% legend('tautest5','taulist5')
% hold on;
% subplot(2,3,6)
% plot(1:n,tau_test(:,6),'b',1:n,tau_real(:,6),'r');
% xlabel('time（s）');
% ylabel('torque（N/m）');
% legend('tautest6','taulist6')
% % xlswrite('tau6',tau_test(:,6),'A');
% % xlswrite('tau6',tau_real(:,6),'B');
% frame = getframe(A);
%% 误差
% figure(2)
% subplot(2,3,1)
% plot(1:n,tau_test(:,1)-tau_real(:,1));
% legend('tau_test1-tau_real1');
% hold on;
% subplot(2,3,2)
% plot(1:n,tau_test(:,2)-tau_real(:,2));
% legend('tau_test-tau_real2');
% hold on;
% subplot(2,3,3)
% plot(1:n,tau_test(:,3)-tau_real(:,3));
% legend('tau_test3-tau_real3');
% hold on;
% subplot(2,3,4)
% plot(1:n,tau_test(:,4)-tau_real(:,4));
% legend('tau_test4-tau_real4');
% hold on;
% subplot(2,3,5)
% plot(1:n,tau_test(:,5)-tau_real(:,5));
% legend('tau_test5-tau_real5');
% hold on;
% subplot(2,3,6)
% plot(1:n,tau_test(:,6)-tau_real(:,6));
% legend('tau_test6-tau_real6');
end