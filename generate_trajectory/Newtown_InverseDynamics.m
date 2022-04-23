function [taulist]= Newtown_InverseDynamics(thetalist, dthetalist, ddthetalist, g,...
                                   dh_list, mass_list, mass_center_list, inertia_tensor_list, f_tip)
%---------��������-----------------
% thetalist:6x1,�ؽڱ�����dthetalist:6x1,�ؽڱ���һ�׵���, ddthetalist:6x1,�ؽڱ������׵���
% g: 1x1,�������ٶ�
% dh_list:6x4,modified_DH����
% mass_list:6x1,����������mass_center_list:6x3,�����������������ϵ{i}���꣬
% inertia_tensor_list: 3x3x6,�����������������ϵ�Ĺ�������
% f_tip: 2x3,��е��ĩ��ʩ������������

%taulist:6x1,���ؽ���������


%R��3x3x6,��ת����P��3x6,��һ��������ϵ��ǰһ��������ϵ�е�λ��
%w��3x6,���˽��ٶȣ�dw��3x6,���˽Ǽ��ٶȣ�dv��3x6,�����߼��ٶȣ�dvc��3x6,���������߼��ٶ�
%Ic:3x3x6,��ͬ��inertia_tensor_list
%Pc:3x6, mass_center_list��ת��
%F:3x6,�������ܺ�����N:3x6,�������ܺ�����
%f:3x6,ǰһ�����һ�������n:3x6,ǰһ�����һ�������

dof_num = size(dthetalist,1);

alpha = dh_list(:,1);
a = dh_list(:,2);
d = dh_list(:,3);
theta = dh_list(:,4);

m = mass_list;
Pc = mass_center_list';
Ic = inertia_tensor_list;

Z=[0;0;1];
%ת��������

theta = theta + thetalist;
T=zeros(4,4,dof_num);R=zeros(3,3,dof_num);P=zeros(3,dof_num);
for i=1:dof_num
    T(:,:,i)=[cos(theta(i))                   -sin(theta(i))                0                  a(i)
              sin(theta(i))*cos(alpha(i))  cos(theta(i))*cos(alpha(i))  -sin(alpha(i))   -d(i)*sin(alpha(i))
              sin(theta(i))*sin(alpha(i))  cos(theta(i))*sin(alpha(i))  cos(alpha(i))    d(i)*cos(alpha(i))
              0                              0                              0                  1];
    R(:,:,i)=T(1:3,1:3,i);
    P(:,i)=T(1:3,4,i);
end

TT = eye(4,4);
for i = 1:dof_num
    TT = TT*T(:,:,i);
end

%�˶�ѧ�������
w0 = zeros(3,1); dw0 = zeros(3,1);
dv0 = [0;0;g];
w = zeros(3,dof_num); dw = zeros(3,dof_num);
dv = zeros(3,dof_num); dvc = zeros(3,dof_num);
F = zeros(3,dof_num); N = zeros(3,dof_num);

%i = 0
w(:,1) = R(:,:,1)' * w0 + dthetalist(1) * Z;
dw(:,1) = R(:,:,1)' * dw0 + cross(R(:,:,1)' * w0, dthetalist(1) * Z) + ddthetalist(1) * Z;
dv(:,1) = R(:,:,1)' * (cross(dw0,P(:,1)) + cross(w0,cross(w0, P(:,1))) + dv0);
dvc(:,1) = cross(dw(:,1), Pc(:,1))+cross(w(:,1), cross(w(:,1), Pc(:,1))) + dv(:,1);
for i = 1:dof_num-1
   w(:,i+1) = R(:,:,i+1)' * w(:,i) + dthetalist(i+1) * Z ;
   dw(:,i+1) = R(:,:,i+1)' * dw(:,i) + cross(R(:,:,i+1)' * w(:,i), dthetalist(i+1) * Z)+ ddthetalist(i+1) * Z;
   dv(:,i+1) = R(:,:,i+1)' * (cross(dw(:,i), P(:,i+1)) + cross(w(:,i), cross(w(:,i), P(:,i+1))) + dv(:,i));
   dvc(:,i+1) = cross(dw(:,i+1), Pc(:,i+1)) + cross(w(:,i+1), cross(w(:,i+1), Pc(:,i+1))) + dv(:,i+1);
end

for i = 1:dof_num
   F(:,i)=m(i)*dvc(:,i) ;
   N(:,i)=Ic(:,:,i) * dw(:,i) + cross(w(:,i), Ic(:,:,i) * w(:,i));
end

%����ѧ�������
%�ȼ����6����������
taulist = zeros(dof_num,1);
f=zeros(3,dof_num); n=zeros(3,dof_num);

f(:,dof_num) = F(:,dof_num) + f_tip(1,:)';
n(:,dof_num) = N(:,dof_num) + f_tip(2,:)' + cross(Pc(:,dof_num), F(:,dof_num));
taulist(dof_num) = n(:,dof_num)' * Z;
%�ټ����5��1����������
for i=dof_num-1:-1:1
   f(:,i) = R(:,:,i+1) * f(:,i+1) + F(:,i);
   n(:,i) = N(:,i) + R(:,:,i+1) * n(:,i+1) + cross(Pc(:,i), F(:,i))...
            + cross(P(:,i+1), R(:,:,i+1) * f(:,i+1));
   taulist(i) = n(:,i)' * Z;
   
end     

end