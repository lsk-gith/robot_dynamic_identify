function [Bp,PS] = calculteBp(op1,op2,byta)
%%%求最小惯性参数集的sym解
%%%%%%%%%质心位置%%%%%%%%%%%%%%%%%%%%%%%
Pc = sym('Pc',[3,6]);
assume(Pc,'real');
Rv = sym('Rv',[6,1]);
assume(Pc,'real');
Rc = sym('Rc',[6,1]);
assume(Pc,'real');
Rf = sym('Rf',[6,1]);
assume(Pc,'real');
%%%%%%%%%%%%质量%%%%%%%%%%%%%%%%%%%%%%%
m = sym('m',[6,1]);
assume(m, 'real');
%%%%%%%%%%%%转动惯量%%%%%%%%%%%%%%%%%%%%
Ic(:,:,1) = sym('Ic1',[3,3]);
Ic(:,:,2) = sym('Ic2',[3,3]);
Ic(:,:,3) = sym('Ic3',[3,3]);
Ic(:,:,4) = sym('Ic4',[3,3]);
Ic(:,:,5) = sym('Ic5',[3,3]);
Ic(:,:,6) = sym('Ic6',[3,3]);
assume(Ic, 'real');
for i=1:6
    Ic(1,2,i)=Ic(2,1,i);
    Ic(1,3,i)=Ic(3,1,i);
    Ic(2,3,i)=Ic(3,2,i);
end
P = [];
for i=1:6
 P = [P,m(i),m(i)*Pc(1,i),m(i)*Pc(2,i),m(i)*Pc(3,i),Ic(1,1,i),Ic(1,2,i),Ic(1,3,i),Ic(2,2,i),Ic(2,3,i),Ic(3,3,i)];
end
F = [];
for i=1:6
  F = [F,Rv(i),Rc(i)];
%   F = [F,Rv(i),Rc(i),Rf(i)];
end
 PS = [P,F]';
 P = P';
 kk2 = [];
 for i = 1:size(op1,2)
    kk2 = [kk2,P(op1(i))];
 end
kk = [];
for i = 1:size(op2,2)
    kk = [kk,P(op2(i))];
end
Bp = kk' + byta*kk2';
Bp = [Bp;F'];
end