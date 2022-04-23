function [A] = Compute_Amatrix(R,P)
%COMPUTE_AMATIX 
%nΪ�ؽ���
%R��3x3xn,��ת����P��3xn,��һ��������ϵ��ǰһ��������ϵ�е�λ��
%A��6x6xn,
n = size(R,3);
A = zeros(6,6,n);
Px = zeros(3,3,n);

for i = 1:n
    Px(:,:,i)=[0 -P(3,i) P(2,i);P(3,i) 0 -P(1,i);-P(2,i) P(1,i) 0]; 
end

A(:,:,n) = eye(6);
for i = 1:n-1
    A(:,:,i) = [R(:,:,i+1) zeros(3,3); Px(:,:,i+1)*R(:,:,i+1) R(:,:,i+1)];
end 

end
