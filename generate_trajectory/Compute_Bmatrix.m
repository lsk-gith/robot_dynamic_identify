function [B] = Compute_Bmatrix(w,dw,dv)
%COMPUTE_BMATIX
%n为关节数
%w：3xn,连杆角速度，dw：3xn,连杆角加速度，dv：3xn,连杆线加速度，
%B：6x10xn,
n = size(w,2);
B = zeros(6,10,n);
wx = zeros(3,3,n); dwx = zeros(3,3,n); dvx = zeros(3,3,n);
w_dot = zeros(3,6,n); dw_dot = zeros(3,6,n);

for i = 1:n
    wx(:,:,i) = [0 -w(3,i) w(2,i); w(3,i) 0 -w(1,i); -w(2,i) w(1,i) 0]; 
    dwx(:,:,i) = [0 -dw(3,i) dw(2,i); dw(3,i) 0 -dw(1,i); -dw(2,i) dw(1,i) 0];
    dvx(:,:,i) = [0 -dv(3,i) dv(2,i); dv(3,i) 0 -dv(1,i); -dv(2,i) dv(1,i) 0];
    
    w_dot(:,:,i) = [w(1,i) w(2,i) w(3,i) 0 0 0;0 w(1,i) 0 w(2,i) w(3,i) 0;0 0 w(1,i) 0 w(2,i) w(3,i)];
    dw_dot(:,:,i) = [dw(1,i) dw(2,i) dw(3,i) 0 0 0;0 dw(1,i) 0 dw(2,i) dw(3,i) 0;0 0 dw(1,i) 0 dw(2,i) dw(3,i)];
end

for i = 1:n
    B(:,:,i) = [dv(:,i) dwx(:,:,i)+wx(:,:,i)*wx(:,:,i) zeros(3,6);zeros(3,1) -dvx(:,:,i) dw_dot(:,:,i)+wx(:,:,i)*w_dot(:,:,i)];
end 

end
