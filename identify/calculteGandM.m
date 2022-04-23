function [G,M]=calculteGandM(dh,BP,op2,P)
q = sym('q',[6,1]);
assume(q,'real');
% dq = sym('dq',[6,1]);
% assume(dq,'real');
% ddq = sym('ddq',[6,1]);
% assume(ddq,'real');
dq=[0,0,0,0,0,0]';
ddq=[0,0,0,0,0,0]';
Bp = sym('Bp',[size(BP,1)-12,1]);
assume(Bp,'real');
K = calculteBigK(dh,q,dq,ddq);
kk = [];
for i = 1:size(op2,2)
    kk = [kk,K(:,op2(i))];
end
G = kk * P * Bp;%6X1
for i=1:6
    ddq(i)=1;
    K = calculteBigK(dh,q,dq,ddq);
    kk = [];
    for ii = 1:size(op2,2)
        kk = [kk,K(:,op2(ii))];
    end
    MG=  kk * P * Bp;%6X1
    MG = MG - G;
    for j=1:6
        M(j,i)=MG(j);
        ddq(i)=0;
    end
end
end