function [w,dw,dv] = newton_euler_dynamic(dh,theta,dtheta,ddtheta)
%�κ����Ǽ���������˵��˶�����������ʱ�Ľ�dh�������ؽڽǶȲ�����theta,dtheta,ddtheta��������ϵ֮��Ĺ�ϵ��ת����ϵ��R����ƽ����ϵ��P�������Ϊ�������˵��˶��������������ٶȣ�w�����Ǽ��ٶȣ�dw�����ͼ��ٶȣ�dv��
number = size(dh,1);
[R ,P]= compute_frame_transform(dh,theta);
g = 9.81;
Z = [0,0,1]';
w0 = zeros(3,1); dw0 = zeros(3,1);dv0 = [0;0;g];
% w = zeros(3,number); 
% dw = zeros(3,number);
% dv = zeros(3,number);
%% 1-n���ƹ�ʽ
%��һ�ؽ�
w(:,1) = R(:,:,1)' * w0 + dtheta(1) * Z;
dw(:,1) = R(:,:,1)' * dw0 + cross(R(:,:,1)' * w0, dtheta(1) * Z) + ddtheta(1) * Z;
dv(:,1) = R(:,:,1)' * (cross(dw0,P(:,1)) + cross(w0,cross(w0, P(:,1))) + dv0);
%����n-1�ؽ�
for i = 1:number-1
   w(:,i+1) = R(:,:,i+1)' * w(:,i) + dtheta(i+1) * Z ;
   dw(:,i+1) = R(:,:,i+1)' * dw(:,i) + cross(R(:,:,i+1)' * w(:,i), dtheta(i+1) * Z)+ ddtheta(i+1) * Z;
   dv(:,i+1) = R(:,:,i+1)' * (cross(dw(:,i), P(:,i+1)) + cross(w(:,i), cross(w(:,i), P(:,i+1))) + dv(:,i));
end

end
