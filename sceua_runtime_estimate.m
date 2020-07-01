% Estimating the computational requirements for running the SCE-UA
% optimization algorithm
%
% 11/20/2019

N = @(n,p) 2*n.^2+2.*n.*p+3.*n+p+2;
p = 0:12;
y_constant_n = N(4, p);

n = 0:12;
y_constant_p = N(n, 3);

%%
figure

subplot(1,2,1)
plot(n,y_constant_p, 'linewidth', 3)
grid on
title('N_{complexes} = 4')
% ylim([0,400])
xlabel('Number of parameters')
ylabel('Function calls per iteration')
set(gca, 'fontsize', 32)

subplot(1,2,2)
plot(p,y_constant_n, 'linewidth', 3)
grid on
title('N_{pars} = 3')
% ylim([0,400])
xlabel('Number of complexes')
ylabel('Function calls per iteration')
set(gca, 'fontsize', 32)
