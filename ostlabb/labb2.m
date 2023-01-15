clear all
close all
clc

% particles and positions
dim = 2;
NP = 2;
m = [1,1];
M = diag(m);
L = 1;
x = [0,0 ; 1.8,0];

% springs
NS = 1;
ks = 10;
kd = 0.5;

% velocity
vini = 0.01;
v = [0,5;0,-5];

%timestep - environment variables
g = 0;
timesteps = 1500;
dt = 0.01;
time = 0:dt:timesteps*dt;
N = length(time);

%grafik
daspect([1,1,1]);   % Skalar axlarna lika
axis([-5,5,-5,5]);    % Axlarnas intervall
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
    spring(spring_number).KS=ks; % spring coefficient
    spring(spring_number).KD=kd; % damping coefficient

end
for i=1:NS
    LINJE_a1(i)=line([x(spring(i).from,1),x(spring(i).to,1)],[x(spring(i).from,2),x(spring(i).to,2)]);
end

Angularmomentum = zeros(1,N);
Springlength = zeros(1,N);
Angularfreq = zeros(1,N);
Es = zeros(1,N);
Ek = zeros(1,N);
Ep = zeros(1,N);

F = zeros(NP,dim)

%%

for k = 1:NS

    r = x(spring(k).from,:)-x(spring(k).to,:)
    rdot = v(spring(k).from,:)-v(spring(k).to,:);

    f = -((spring(k).KS)*(norm(r)-spring(k).length)+spring(k).KD*(dot(rdot,r)/norm(r))).*r/norm(r);

    F(spring(k).from,:) = F(spring(k).from,:) + f;
    F(spring(k).to,:) = F(spring(k).to,:) - f;

    Es = Es + spring(k).KS*(norm(r)-spring(k).length).^2./2;
end

v = v - dt*inv(M)*F;%./2;

%%

for n = 1:N

    F = zeros(NP,dim)
    for k = 1:NS

        r = x(spring(k).from,:)-x(spring(k).to,:)
        rdot = v(spring(k).from,:)-v(spring(k).to,:);

        f = -((spring(k).KS)*(norm(r)-spring(k).length)+spring(k).KD*(dot(rdot,r)/norm(r))).*r/norm(r);

        F(spring(k).from,:) = F(spring(k).from,:) + f;
        F(spring(k).to,:) = F(spring(k).to,:) - f;

        Es = Es + spring(k).KS*(norm(r)-spring(k).length).^2./2;
    end
    F(:,dim) = F(:,dim)-m'*g;
    v_prev = v;
    v = v + dt*inv(M)*F;
    x = x + dt*v;
    Angularmomentum(n) = m*(x(:,1).*v(:,2)-x(:,2).*v(:,1))
    %Angularfreq(n) = norm(v(1,:)./(norm(r)./2))
    Ek(n) = sum(m.*(sum(((v+v_prev)./2).^2)./2))
    Ep(n) = sum(g.*m.*x(:,2)')

    %uppdaterad grafik
    for i = 1:NS
        set(LINJE_a1(i),'XData',[x(spring(i).from,1),x(spring(i).to,1)],'YData',[x(spring(i).from,2),x(spring(i).to,2)]);
    end
    for i=1:NP
        set(BOLL_a(i),'Position',[x(i,1)-Ra,x(i,2)-Ra,2*Ra,2*Ra]);
    end

    pause(dt/1000);
end

Etot = Es + Ek + Ep

figure(2)
plot(Etot)