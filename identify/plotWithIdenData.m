load('tauAndTauIden.mat');
Fram = plotverify(tau,tauIden);
S = calculateS(tau,tauIden);
im1=frame2im(Fram);