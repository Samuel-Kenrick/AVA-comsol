function find_max_FMM
% This program will find the maximum amplitude of u for both omega=pi
% and omega= 1.5*pi from data derived via CMP in Computerlab 2 in
% the course "Mathematics of Physical Models B, 10.5 ECTS" given at
% UmU.

%% Loading data
% u omega#(:,1) defines the geometry, u omega1(:,2:end) is the
% computed data.
u_omega1 = load('omega1.txt');
u_omega2 = load('omega2.txt');


 %% Finding maximum
 s = size(u_omega1);
 max_omega1(1) = 0;
 max_omega2(1) = 0;
 Max_omega1 = max_omega1(1);
 Max_omega2 = max_omega2(1);
 for n = 2 : s(2)
 max_omega1(n) = max(abs(u_omega1(:,n)));
 max_omega2(n) = max(abs(u_omega2(:,n)));

 if max_omega1(n)>Max_omega1
 Max_omega1 = max_omega1(n);
 MAX_n1 = n;
 end

 if max_omega2(n)>Max_omega2
 Max_omega2 = max_omega2(n);
 MAX_n2 = n;
 end
 end

 %% Ploting the result
 plot(u_omega1(:,1),u_omega1(:,MAX_n1),u_omega2(:,1),u_omega2(:,MAX_n2))
 xlabel('x−coordinate [m]');
 ylabel('Amplitude [m]');
 legend(['\omega=\pi, ', sprintf('t=%.1f',0.1*(MAX_n1-2))],['\omega=1.5\pi, ',...
sprintf('t=%.1f',0.1*(MAX_n2-2))],'fontsize',11)