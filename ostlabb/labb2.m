clear all
close all
clc

% particles and positions
NP = 2;
m = [1,1];
M = diag(m);
L = 1;
x = [0,0 ; 1.8,0];

% springs
NS = 1;
ks = 10;
kd = 0;

% velocity
vini = 0.01;
v = [0;0];

%timestep
timesteps = 1500;
dt = 0.1;
time = 0:dt:timesteps*dt;
N = length(time);

%grafik
daspect([1,1,1]);   % Skalar axlarna lika
axis([-3,3,-3,3]);    % Axlarnas intervall
hold on
Ra = 0.05;
for i=1:NP
    BOLL_a(i)=rectangle('Position',[x(i,1)-Ra,x(i,2)-Ra,2*Ra,2*Ra],'Curvature',[1,1],'EdgeColor','r','FaceColor','w');
end

spring_number=0;
for i=1:NS % loop over springs

spring_number=spring_number+1;
spring(spring_number).from=i; % number of the ''from'' particle
spring(spring_number).to=i+1; % number of the ''to'' particle
spring(spring_number).length=L; % spring rest length
spring(spring_number).KS=KS; % spring coefficient
spring(spring_number).KD=KD; % damping coefficient
end
for i=1:NS
LINJE_a1(i)=line([x(spring(i).from,1),x(spring(i).to,1)],[x(spring(i).from,2),x(spring(i).to,2)]);
end

Ek = zeros(length(m),length(tspan));
Ep = zeros(1,length(tspan));
Es = zeros(1,length(tspan));
angularMomentum = zeros(1,length(tspan));
R = zeros(1,length(tspan));

for n = 1:N
    fspring = -ks(x-L)-kd*v;
    F = 0.5*(0.5-rand(NP,2))-0.1*v;
    v = v + dt*inv(M)*F;
    x = x + dt*v;

    %uppdaterad grafik
    for i=1:NP
        set(BOLL_a(i),'Position',[x(i,1)-Ra,x(i,2)-Ra,2*Ra,2*Ra]);
    end
    set(LINJE1,'XData',[x(1,1),x(2,1)],'YData',[x(1,2),x(2,2)]);
    pause(0.1);
end