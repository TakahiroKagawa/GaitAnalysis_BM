function []=BalanceMap(FigureNum, y_lim, yd_lim, Omega0)
%Å@This function draws a balance map figure
%  

if nargin < 1
    FigureNum=2;
end
if nargin < 2
    load('balancemap_limit.mat');
    Omega0 = par.ChiOmega0;
end


figure(FigureNum);hold on
% phase
psi=-3:0.01:3;

% Contour plots of energy ratio
energy_p=[0.25:0.25:4];
energy_n=[0.25:0.25:4];
for cnt1=1:length(energy_p);
    chi=sqrt(energy_p(cnt1))/Omega0*sinh(psi);
    chid=sqrt(energy_p(cnt1))*cosh(psi);
    figure(FigureNum);
    plot(Omega0*chi,chid,'-k');
    plot(Omega0*chi,-chid,'-k');
end
% Contour plots of phase difference
plot([0 0],[-3 3],'-k');
plot([-3 3],[0 0],'-k');
chi=[-3 3];
plot(Omega0*chi,1/tanh(0.5)*Omega0*chi,'-k');
plot(Omega0*chi,1/tanh(1)*Omega0*chi,'-k');
plot(Omega0*chi,1/tanh(1.5)*Omega0*chi,'-k');

plot(Omega0*chi,-1/tanh(0.5)*Omega0*chi,'-k');
plot(Omega0*chi,-1/tanh(1)*Omega0*chi,'-k');
plot(Omega0*chi,-1/tanh(1.5)*Omega0*chi,'-k');


plot(Omega0*chi,tanh(0.5)*Omega0*chi,'-k');
plot(Omega0*chi,tanh(1)*Omega0*chi,'-k');
plot(Omega0*chi,tanh(1.5)*Omega0*chi,'-k');
plot(Omega0*chi,-tanh(0.5)*Omega0*chi,'-k');
plot(Omega0*chi,-tanh(1)*Omega0*chi,'-k');
plot(Omega0*chi,-tanh(1.5)*Omega0*chi,'-k');

% boundary lines of fall states
for cnt1=1:length(energy_n);
    chi=sqrt(energy_n(cnt1))/Omega0*cosh(psi);
    chid=sqrt(energy_n(cnt1))*sinh(psi);
    figure(FigureNum);
    plot(Omega0*chi,chid,'-k');
    plot(-Omega0*chi,chid,'-k');
end



figure(FigureNum);
plot([-2 2],[-2 2],'-k');
plot([-2 2],[2 -2],'-k','LineWidth',2);
plot(Omega0*y_lim,yd_lim,'-k','Linewidth',2);

axis([-1.0 1.0 0 1.2])
figure(FigureNum);
set(gca,'FontSize',16);
xlabel('Position (m)');
ylabel('Velocity (m/s)')

