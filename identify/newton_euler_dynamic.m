function [w,dw,dv] = newton_euler_dynamic(dh,theta,dtheta,ddtheta)
%次函数是计算各个连杆的运动参数，输入时改进dh参数，关节角度参数（theta,dtheta,ddtheta）和坐标系之间的关系（转动关系：R，和平动关系：P），输出为各个连杆的运动参数，包括角速度（w），角加速度（dw），和加速度（dv）
number = size(dh,1);
[R ,P]= compute_frame_transform(dh,theta);
g = 9.81;
Z = [0,0,1]';
w0 = zeros(3,1); dw0 = zeros(3,1);dv0 = [0;0;g];
% w = zeros(3,number); 
% dw = zeros(3,number);
% dv = zeros(3,number);
%% 1-n外推公式
%第一关节
w(:,1) = R(:,:,1)' * w0 + dtheta(1) * Z;
dw(:,1) = R(:,:,1)' * dw0 + cross(R(:,:,1)' * w0, dtheta(1) * Z) + ddtheta(1) * Z;
dv(:,1) = R(:,:,1)' * (cross(dw0,P(:,1)) + cross(w0,cross(w0, P(:,1))) + dv0);
%后面n-1关节
for i = 1:number-1
   w(:,i+1) = R(:,:,i+1)' * w(:,i) + dtheta(i+1) * Z ;
   dw(:,i+1) = R(:,:,i+1)' * dw(:,i) + cross(R(:,:,i+1)' * w(:,i), dtheta(i+1) * Z)+ ddtheta(i+1) * Z;
   dv(:,i+1) = R(:,:,i+1)' * (cross(dw(:,i), P(:,i+1)) + cross(w(:,i), cross(w(:,i), P(:,i+1))) + dv(:,i));
end

end
