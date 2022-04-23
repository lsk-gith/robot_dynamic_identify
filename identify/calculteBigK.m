function [KK,FF] = calculteBigK(dh,q,dq,ddq)
KK = [];
FF=[];
for i=1:size(q,2)
[R,P] = compute_frame_transform(dh,q(:,i));
[w,dw,dv] = newton_euler_dynamic(dh,q(:,i), dq(:,i), ddq(:,i));
[K,F]= cacluate_K(dh,w,dw,dv,R,P,q(:,i),dq(:,i), ddq(:,i));
FF = [FF;F];
KK = [KK;K];
end
end