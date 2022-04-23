function [KK]=calculteBB(dh,G,P,CC,BP,op2)
%%∑÷¿Îø∆ œ¡¶
q = sym('q',[6,1]);
assume(q,'real');
ddq=[0,0,0,0,0,0]';
dq = [0,0,0,0,0,0]';
Bp = sym('Bp',[size(BP,1)-12,1]);
assume(Bp,'real');
i = 1;
for nn = 1:5
    for ll=(nn+1):6
        dq(nn)=1;dq(ll)=1;
        K = calculteBigK(dh,q,dq,ddq);
        kk = [];
        for ii = 1:size(op2,2)
            kk = [kk,K(:,op2(ii))];
        end
        BCG =  kk * P * Bp;
        B=BCG-G-CC(:,nn)-CC(:,ll);
        KK(:,i)= B;
        dq(nn)=0;dq(ll)=0;
        i = i + 1;
    end
end
end