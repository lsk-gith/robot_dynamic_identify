function [op1,op2,BP,byta,P]=qrLeast(Kk,tau,F)
%%%%%这里需要两次qr分解，第一次是为了选出KK中independent项（R（i,i）=0）, 第二次qr分解是为了得到byta值，求最小惯性参数集里面的sym量(这次qr分解的对象是调整过的KK(independent,dependent))
[~,RV]=qr(Kk);
op1=[];
op2=[];
for i = 1:size(RV,2)
    if(abs(RV(i,i))<(1e-5))
        op1 = [op1,i];%%linerly dependent
    else
        op2 = [op2,i];%%independent
    end
end
kk2 = [];
for i = 1:size(op1,2)
    kk2 = [kk2,Kk(:,op1(i))];
end
kk = [];
for i = 1:size(op2,2)
    kk = [kk,Kk(:,op2(i))];
end

[Q,R]=qr([kk,kk2]);
B = rank(Kk);
R1 = R(1:B,1:B);
R2 = R(1:B,B+1:end);
Q1 = Q(:,1:B);
W1 = Q1 * R1;
norm(W1-kk)%%%%验证W1和kk拆分正确性
byta = R1^-1*R2; %%最小惯性参数集就是kk对应的惯性参数集里面的惯性参数加上byta*kk2对应的惯性参数集里面的惯性参数
% fprintf('####');
% norm(kk*byta-kk2)

%%%方法1
% KK = [Kk(:,10),Kk(:,12:13),Kk(:,15:end)];
% [Q,R,P]=qr(KK);
% % R3 = R(1:B,1:B);
% % fprintf('******');
% % norm(R1-R3)
% KK = KK * P;
% BP = (KK'*KK)^-1*KK'*tau;


%%这里还不能使用最小二乘法做辨识，因为rank(kk)秩和rank([kk,tauFilt((startPoint1-1)*6+1:(endPoint1)*6)])秩不相同，这里要引入另一个量error vector
%%% tau = kk * Bp + row  (剩下的工组就是怎么确定row),这个误差矩阵是电机产生的和电机采集数据有关，
[Q R P] = qr(kk);
norm(kk*P-Q*R)
kk = kk * P;
kk = [kk,F];
rank(kk)
rank([kk,tau])
BP = (kk'*kk)^-1*kk'*tau;

end