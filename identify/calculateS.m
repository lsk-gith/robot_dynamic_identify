function S=calculateS(tau,tauiden)
%%求解总体标准差
%%tau 是没有化成6n * 1的形式
n = size(tau,1);
s = [];
for i =1:n
    s = [s;(tau(i)-tauiden(i))^2];
end

for i = 1:6
    for j = 1:n/6
        ss(j,i) = s(6*(j-1)+i);
    end
end
for i =1:6
    
S(i) = (sum(ss(:,i))/(n/6))^0.5;

end
end