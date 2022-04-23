function [K,F] = cacluate_K(dh,w,dw,dv,R,P,q,dq, ddq)
%%次函数是参数辨识中用于计算最小惯性参数集分离出来的右乘部分（K），输入为改进dh坐标系（dh),各个连杆的运动参数(角速度:w,角加速度:dw,和加速度dv),坐标系之间的关系（转动关系：R，和平动关系：P）
n = size(dh,1);
switch class(q(1))
    case 'double'
        RI1 = [];DK = [];
    otherwise
end

for i = 1:n
    px(:,:,i) = [0,-P(3,i),P(2,i);P(3,i),0,-P(1,i);-P(2,i),P(1,i),0];
    wx(:,:,i) = [0,-w(3,i),w(2,i);w(3,i),0,-w(1,i);-w(2,i),w(1,i),0];
    dwx(:,:,i) = [0,-dw(3,i),dw(2,i);dw(3,i),0,-dw(1,i);-dw(2,i),dw(1,i),0];
    dvx(:,:,i)=[0,-dv(3,i),dv(2,i);dv(3,i),0,-dv(1,i);-dv(2,i),dv(1,i),0];
    Wdot(:,:,i) = [w(:,i)',zeros(1,3);0,w(1,i),0,w(2,i),w(3,i),0;0,0,w(1,i),0,w(2,i),w(3,i)];
    dWdot(:,:,i) = [dw(:,i)',zeros(1,3);0,dw(1,i),0,dw(2,i),dw(3,i),0;0,0,dw(1,i),0,dw(2,i),dw(3,i)];
    DK(:,:,i) = [dv(:,i),dwx(:,:,i)+wx(:,:,i)*wx(:,:,i),zeros(3,6); ...
        zeros(3,1),-dvx(:,:,i),dWdot(:,:,i)+wx(:,:,i)*Wdot(:,:,i)];
end
for i=1:n-1
     RI1(:,:,i) = [R(:,:,i+1),zeros(3,3);px(:,:,i+1)*R(:,:,i+1),R(:,:,i+1)];
end
%%  计算大K
DKKKK = [];
for i =1:n-1
    D = [];
    RR = eye(6,6);
    for m = i:n-1
        RR = RR * RI1(:,:,m);
        D = [D,RR * DK(:,:,m+1)];
    end
    DKKK = [zeros(6,(i-1)*10),DK(:,:,i),D];
    DKKKK = [DKKKK;DKKK];
end
DKKKK = [DKKKK;zeros(6,(n-1)*10),DK(:,:,n)];
% size(DKKKK)

%取tau值
K = [];
F = [];
for i =1:n
    K = [K;DKKKK(6*i,:)];
end
% size(K)
%% 加摩擦力

for i=1:n
 %    F = [F;zeros(1,2*(i-1)),dq(i),sign(dq(i)),zeros(1,10-2*(i-1))];
    F = [F;zeros(1,2*(i-1)),dq(i),tanh(dq(i)/0.001),zeros(1,10-2*(i-1))];
%     F = [F;zeros(1,3*(i-1)),dq(i),sign(dq(i)),ddq(i),zeros(1,15-3*(i-1))];
%     F = [F;zeros(1,3*(i-1)),dq(i),tanh(dq(i)/0.05),ddq(i),zeros(1,15-3*(i-1))];

%     F = [F;zeros(1,3*(i-1)),1,1,1,zeros(1,15-3*(i-1))];
end
% K = [K,F];


end