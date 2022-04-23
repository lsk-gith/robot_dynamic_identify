function [J,K,Kf] = Compute_Dynmatrix(thetalist, dthetalist, ddthetalist, g, dh_list, nf)
% 动力学模型线性化,计算线性回归矩阵
% 输入: thetalist: nx1,关节角度; dthetalist: nx1,关节角速度; dthetalist: nx1,关节角加速度; 
%       g: 1x1,重力加速度; dh_list: nx4, 修正DH参数; nf: 1x1, 摩擦力项数
% 输出: J: nx6, 与外力有关的线性回归矩阵
%       K: nx(10*n), 与机器人本体参数有关的线性回归矩阵
%       Kf: nx(nf*n), 与摩擦力有关的线性回归矩阵,
%           nf=2: 库伦粘滞摩擦模型;  nf=4: 正反库伦粘滞摩擦模型; nf=5: 正反库伦粘滞摩擦模型+电机惯性补偿; 
%%
n = size(dh_list,1);

alpha = dh_list(:,1);
a = dh_list(:,2);
d = dh_list(:,3);
theta = dh_list(:,4);

Z=[0;0;1];
%转换矩阵建立

theta = theta + thetalist;
T=zeros(4,4,n);R=zeros(3,3,n);P=zeros(3,n);
for i=1:n
    T(:,:,i)=[cos(theta(i))                   -sin(theta(i))                0                  a(i)
              sin(theta(i))*cos(alpha(i))  cos(theta(i))*cos(alpha(i))  -sin(alpha(i))   -d(i)*sin(alpha(i))
              sin(theta(i))*sin(alpha(i))  cos(theta(i))*sin(alpha(i))  cos(alpha(i))    d(i)*cos(alpha(i))
              0                              0                              0                  1];
    R(:,:,i)=T(1:3,1:3,i);
    P(:,i)=T(1:3,4,i);
end

%运动学正向递推
w0 = zeros(3,1); dw0 = zeros(3,1);
dv0 = [0;0;g];
w = zeros(3,n); dw = zeros(3,n);
dv = zeros(3,n);

%i = 0
w(:,1) = R(:,:,1)' * w0 + dthetalist(1) * Z;
dw(:,1) = R(:,:,1)' * dw0 + cross(R(:,:,1)' * w0, dthetalist(1) * Z) + ddthetalist(1) * Z;
dv(:,1) = R(:,:,1)' * (cross(dw0,P(:,1)) + cross(w0,cross(w0, P(:,1))) + dv0);
for i = 1:n-1
   w(:,i+1) = R(:,:,i+1)' * w(:,i) + dthetalist(i+1) * Z ;
   dw(:,i+1) = R(:,:,i+1)' * dw(:,i) + cross(R(:,:,i+1)' * w(:,i), dthetalist(i+1) * Z)+ ddthetalist(i+1) * Z;
   dv(:,i+1) = R(:,:,i+1)' * (cross(dw(:,i), P(:,i+1)) + cross(w(:,i), cross(w(:,i), P(:,i+1))) + dv(:,i));
end

A = Compute_Amatrix(R,P);
% A1 = A(:,:,1); A2 = A(:,:,2); A3 = A(:,:,3); A4 = A(:,:,4); A5 = A(:,:,5); A6 = A(:,:,6);
B = Compute_Bmatrix(w,dw,dv);
% B1 = B(:,:,1); B2 = B(:,:,2); B3 = B(:,:,3); B4 = B(:,:,4); B5 = B(:,:,5); B6 = B(:,:,6);


%计算J
V = []; 
for i = 1:n
    V_temp = eye(6,6);
    for j = i:n
    V_temp = V_temp * A(:,:,j);
    end
    V = [V;V_temp];
end

J = [];
for i = 1:n
    J = [J;V(6*i,:)];
end
%计算K
U = [];
for i =1:n-1
    RR = eye(6,6);
    D = [];
    for j = i:n-1
        RR = RR * A(:,:,j);
        D = [D,RR * B(:,:,j+1)];
    end
    U_temp = [zeros(6,(i-1)*10),B(:,:,i),D];
    U = [U;U_temp];
end
U = [U;zeros(6,(n-1)*10),B(:,:,n)]; %最后一次循环要拿出来写

K = [];
for i =1:n
    K = [K;U(6*i,:)];
end

% U11 = B1; U12 = A1*B2; U13 = A1*A2*B3; U14 = A1*A2*A3*B4; U15 = A1*A2*A3*A4*B5; U16 = A1*A2*A3*A4*A5*B6;
% U22 = B2; U23 = A2*B3; U24 = A2*A3*B4; U25 = A2*A3*A4*B5; U26 = A2*A3*A4*A5*B6;
% U33 = B3; U34 = A3*B4; U35 = A3*A4*B5; U36 = A3*A4*A5*B6;
% U44 = B4; U45 = A4*B5; U46 = A4*A5*B6;
% U55 = B5; U56 = A5*B6;
% U66 = B6;
% 
% V = [A1*A2*A3*A4*A5*A6;
%      A2*A3*A4*A5*A6;
%      A3*A4*A5*A6;
%      A4*A5*A6;
%      A5*A6;
%      A6];
%  
% U = [U11 U12 U13 U14 U15 U16;
%     zeros(6,10) U22 U23 U24 U25 U26;
%     zeros(6,10) zeros(6,10) U33 U34 U35 U36;
%     zeros(6,10) zeros(6,10) zeros(6,10) U44 U45 U46;
%     zeros(6,10) zeros(6,10) zeros(6,10) zeros(6,10) U55 U56;
%     zeros(6,10) zeros(6,10) zeros(6,10) zeros(6,10) zeros(6,10) U66];
% 
% % temp = [0 0 0 0 0 1];
% % KK = [temp zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6);
% %      zeros(1,6) temp zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6);
% %      zeros(1,6) zeros(1,6) temp zeros(1,6) zeros(1,6) zeros(1,6);
% %      zeros(1,6) zeros(1,6) zeros(1,6) temp zeros(1,6) zeros(1,6)
% %      zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6) temp zeros(1,6)
% %      zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6) temp] * U;
% %  
% % JJ = [temp zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6);
% %      zeros(1,6) temp zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6);
% %      zeros(1,6) zeros(1,6) temp zeros(1,6) zeros(1,6) zeros(1,6);
% %      zeros(1,6) zeros(1,6) zeros(1,6) temp zeros(1,6) zeros(1,6)
% %      zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6) temp zeros(1,6)
% %      zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6) zeros(1,6) temp] * V; 
%  
% 
% K = [U11(end,:) U12(end,:) U13(end,:) U14(end,:) U15(end,:) U16(end,:);
%     zeros(1,10) U22(end,:) U23(end,:) U24(end,:) U25(end,:) U26(end,:);
%     zeros(1,10) zeros(1,10) U33(end,:) U34(end,:) U35(end,:) U36(end,:);
%     zeros(1,10) zeros(1,10) zeros(1,10) U44(end,:) U45(end,:) U46(end,:);
%     zeros(1,10) zeros(1,10) zeros(1,10) zeros(1,10) U55(end,:) U56(end,:);
%     zeros(1,10) zeros(1,10) zeros(1,10) zeros(1,10) zeros(1,10) U66(end,:)];
% 
% J = [V(6,:);V(12,:);V(18,:);V(24,:);V(30,:);V(36,:)];
%% compute Kf
Kf = zeros(n,nf*n);
epsilon = 0.01;
%Kf = [diag(dthetalist),diag(sign(dthetalist))];
%Kf = [diag(ddthetalist),diag(dthetalist),diag(sign(dthetalist))];
for i = 1:n
    %Kf(i,:) = [zeros(1,nf*(i-1)) [sign(dthetalist(i)) dthetalist(i)] zeros(1,nf*(n-i))];
    if(nf == 2)
        Kf(i,:) = [zeros(1,nf*(i-1)) [tanh(dthetalist(i)/epsilon) dthetalist(i)] zeros(1,nf*(n-i))];
    elseif(nf == 3)
        Kf(i,:) = [zeros(1,nf*(i-1)) [tanh(dthetalist(i)/epsilon) dthetalist(i) ddthetalist(i)] zeros(1,nf*(n-i))];
    elseif(nf == 4)
        Kf(i,:) = [zeros(1,nf*(i-1)) [tanh(dthetalist(i)/epsilon)*(tanh(dthetalist(i)/epsilon)+1)/2  dthetalist(i)*(tanh(dthetalist(i)/epsilon)+1)/2 tanh(dthetalist(i)/0.01)*(1-tanh(dthetalist(i)/epsilon))/2  dthetalist(i)*(1-tanh(dthetalist(i)/epsilon))/2 ] zeros(1,nf*(n-i))];
    elseif(nf == 5)
        Kf(i,:) = [zeros(1,nf*(i-1)) [tanh(dthetalist(i)/epsilon)*(tanh(dthetalist(i)/epsilon)+1)/2  dthetalist(i)*(tanh(dthetalist(i)/epsilon)+1)/2 tanh(dthetalist(i)/0.01)*(1-tanh(dthetalist(i)/epsilon))/2  dthetalist(i)*(1-tanh(dthetalist(i)/epsilon))/2 ddthetalist(i)] zeros(1,nf*(n-i))];
    end
end
end
