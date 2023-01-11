clear all
clc

% particles and positions
NP = 3;
m = [1,1,1];
M = diag(m);
L = 1;
x = [-L/2,0 ; L/2,0 ; 0,L/2];

% velocity
vini = 0.01;
v = [vini,0 ; 0,-vini ; -vini/2,-vini/2];

%timestep
timesteps = 1500;
dt = 0.1;
time = 1:dt:timesteps*dt;
N = length(time);

% m = 1;

%grafik
daspect([1,1,1]);   % Skalar axlarna lika
axis([-3,3,-3,3]);    % Axlarnas intervall
hold on
Ra = 0.1;
for i=1:NP
    BOLL_a(i)=rectangle('Position',[x(i,1)-Ra,x(i,2)-Ra,2*Ra,2*Ra],'Curvature',[1,1],'EdgeColor','r','FaceColor','w');
end
LINJE1=line([x(1,1),x(2,1)],[x(1,2),x(2,2)]);
LINJE2=line([x(3,1),x(2,1)],[x(3,2),x(2,2)]);
LINJE3=line([x(1,1),x(3,1)],[x(1,2),x(3,2)]);

for n = 1:N

    F = 0.5*(0.5-rand(NP,2))-0.1*v;
    v = v + dt*inv(M)*F;
    x = x + dt*v;

    %uppdaterad grafik
    for i=1:NP
        set(BOLL_a(i),'Position',[x(i,1)-Ra,x(i,2)-Ra,2*Ra,2*Ra]);
    end
    set(LINJE1,'XData',[x(1,1),x(2,1)],'YData',[x(1,2),x(2,2)]);
    set(LINJE2,'XData',[x(3,1),x(2,1)],'YData',[x(3,2),x(2,2)]);
    set(LINJE3,'XData',[x(1,1),x(3,1)],'YData',[x(1,2),x(3,2)]);
    pause(0.1); 
end