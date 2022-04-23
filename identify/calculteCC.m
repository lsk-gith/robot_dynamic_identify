function [CC]=calculteCC(dh,G,op2,P,BP)
%%分离离心力项
q = sym('q',[6,1]);
assume(q,'real');
ddq=[0,0,0,0,0,0]';
dq=[0,0,0,0,0,0]';
Bp = sym('Bp',[size(BP,1)-12,1]);
assume(Bp,'real');
for i=1:6
    dq(i)=1;
    K = calculteBigK(dh,q,dq,ddq);
    kk = [];
    for ii = 1:size(op2,2)
        kk = [kk,K(:,op2(ii))];
    end
    CG=  kk * P * Bp;
    CG = CG - G;
    for j=1:6
       CC(j,i)=CG(j);
       dq(i)=0;
    end
end
end