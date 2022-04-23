function [op1,op2,BP,byta,P]=qrLeast(Kk,tau,F)
%%%%%������Ҫ����qr�ֽ⣬��һ����Ϊ��ѡ��KK��independent�R��i,i��=0��, �ڶ���qr�ֽ���Ϊ�˵õ�bytaֵ������С���Բ����������sym��(���qr�ֽ�Ķ����ǵ�������KK(independent,dependent))
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
norm(W1-kk)%%%%��֤W1��kk�����ȷ��
byta = R1^-1*R2; %%��С���Բ���������kk��Ӧ�Ĺ��Բ���������Ĺ��Բ�������byta*kk2��Ӧ�Ĺ��Բ���������Ĺ��Բ���
% fprintf('####');
% norm(kk*byta-kk2)

%%%����1
% KK = [Kk(:,10),Kk(:,12:13),Kk(:,15:end)];
% [Q,R,P]=qr(KK);
% % R3 = R(1:B,1:B);
% % fprintf('******');
% % norm(R1-R3)
% KK = KK * P;
% BP = (KK'*KK)^-1*KK'*tau;


%%���ﻹ����ʹ����С���˷�����ʶ����Ϊrank(kk)�Ⱥ�rank([kk,tauFilt((startPoint1-1)*6+1:(endPoint1)*6)])�Ȳ���ͬ������Ҫ������һ����error vector
%%% tau = kk * Bp + row  (ʣ�µĹ��������ôȷ��row),����������ǵ�������ĺ͵���ɼ������йأ�
[Q R P] = qr(kk);
norm(kk*P-Q*R)
kk = kk * P;
kk = [kk,F];
rank(kk)
rank([kk,tau])
BP = (kk'*kk)^-1*kk'*tau;

end