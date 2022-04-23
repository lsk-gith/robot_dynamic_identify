function [R,P] = compute_frame_transform(dh,thetalist)
%%次函数是用于计算坐标系之间的转换关系，输入为改进dh参数（dh）和角度变化量（thetalist），输出为坐标系之间的转动关系（R），和平动关系（P）。
n = size(dh,1);
alpha = dh(:,1);
a = dh(:,2);
d = dh(:,3);
theta = dh(:,4);
theta = theta + thetalist;
switch class(thetalist(1))
    case 'double'
        T=zeros(4,4,n);R=zeros(3,3,n);P=zeros(3,n);
    otherwise
end

for i=1:n
    T(:,:,i)=[cos(theta(i))                   -sin(theta(i))                0                  a(i)
              sin(theta(i))*cos(alpha(i))  cos(theta(i))*cos(alpha(i))  -sin(alpha(i))   -d(i)*sin(alpha(i))
              sin(theta(i))*sin(alpha(i))  cos(theta(i))*sin(alpha(i))  cos(alpha(i))    d(i)*cos(alpha(i))
              0                              0                              0                  1];
    R(:,:,i)=T(1:3,1:3,i);
    P(:,i)=T(1:3,4,i);
end
end
