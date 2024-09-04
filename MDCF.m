%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------
% Function: calculate the Crossover Frequency by Modal-density method
% Input1: room dimensitions: lx,ly,lz
% Input2: boundary absorption coefficient: alpha 
% Results: MDCF
% --------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_x = 3.2;L_y=1.8*1;L_z=1.4*1;
alpha_values = 0.01:0.001:0.5;

V = L_x*L_y*L_z; % 房间的体积，单位为m^3
S = 2*(L_x*L_y+L_x*L_z+L_z*L_y); % 房间的表面积，单位为m^2
L = 4*(L_x+L_y+L_z);

rho0 = 1.2;
c0 = 343;
zn = rho*c0.*(1+sqrt(1-alpha_values))./(1-sqrt(1-alpha_values));

tau = zn .* V/(2*rho0*c0^2*S);
delta = 1./(2.*tau);

A = (4*pi*V)/(c0^3);
B = (pi*S)/(2*c0^2);
C = (L/(8*c0)) - (3*pi./delta);
discriminant = B^2 - 4*A.*C;
F_CF1 = 1.2*(-B + sqrt(discriminant)) ./ (2*A);

plot(alpha_values,F_CF1,Color='r',LineWidth=2.5,LineStyle='-');hold on;
set(gca,'fontsize',18,'linewidth',2);
ylim([0 2300]); xlim([0 0.5]); grid on;
legend('MDCF-Formula',Location='northwest');