clear all
close all
clc

% environment variables
dim = 2;
Nr = 4;
Nc = 8;
NP = Nr*Nc;
m = zeros(Nr,Nc)+1;
M = diag(m);
L = [1 sqrt(2)];
ks = 100;
kd = 0;
x0 = 1;
y0 = 5;

% velocity
vini = 0.01;
v = [0,5;0,-5];

%timestep - environment variables
g = 1;
timesteps = 1500;
dt = 0.002;
time = 0:dt:timesteps*dt;
N = length(time);

%grafik
daspect([1,1,1]);   % Skalar axlarna lika
axis([0,50,0,50]);    % Axlarnas intervall
hold on
Ra = 0.05;

% Balls
NB = 16;
Ra = 1;


for i=1:NB
    BOLL_b(i)=rectangle('Position',[i*(2*Ra)+0.1*Ra*randn,0,2*Ra,2*Ra],'Curvature',[1,1],'EdgeColor','r','FaceColor','w');
end

%%


for i = 1:Nc*Nr
    x(i,1) = x0 + (Nc-i-1)
    x(i,2) = y0 + (Nr-i-1)
end

%for i = 1:(Nr)
    
%end

for i = 1:NP
    BOLL_a(i)=rectangle('Position',[x0(i,1)-Ra,x0(i,2)-Ra,2*Ra,2*Ra],'Curvature',[1,1],'EdgeColor','#EDB120','FaceColor','#EDB120');       
end

%%
%horisontal
spring_number=0;
particle = 1;
for i=1:(Nr*(Nc-1)) % loop over springs
    
    spring_number=spring_number+1;

    spring(spring_number).from=particle; % number of the ''from'' particle
    spring(spring_number).to=particle+1; % number of the ''to'' particle
    spring(spring_number).length=L(1); % spring rest length
    spring(spring_number).KS=ks; % spring coefficient
    spring(spring_number).KD=kd; % damping coefficient

    particle = particle +1;
    

end

%% verticle
particle = 1
for i=1:(Nr-1)*Nc % loop over springs
    
    spring_number=spring_number+1;

    spring(spring_number).from=particle; % number of the ''from'' particle
    spring(spring_number).to=particle+Nc; % number of the ''to'' particle
    spring(spring_number).length=L(1); % spring rest length
    spring(spring_number).KS=ks; % spring coefficient
    spring(spring_number).KD=kd; % damping coefficient

    particle = particle +1;
  
end

%%
particle = 1
for i=1:(Nr-1)*(Nc-1)+2 % loop over springs to right
    
    spring_number=spring_number+1;

    spring(spring_number).from=particle; % number of the ''from'' particle
    spring(spring_number).to=particle+Nc+1; % number of the ''to'' particle
    spring(spring_number).length=L(2); % spring rest length
    spring(spring_number).KS=ks; % spring coefficient
    spring(spring_number).KD=kd; % damping coefficient

    particle = particle +1;
  
end


%%

particle = 2
for i=1:(Nr-1)*(Nc-1)+(Nr-2) % loop over springs to right
    
    spring_number=spring_number+1;

    spring(spring_number).from=particle; % number of the ''from'' particle
    spring(spring_number).to=particle+Nc-1; % number of the ''to'' particle
    spring(spring_number).length=L(2); % spring rest length
    spring(spring_number).KS=ks; % spring coefficient
    spring(spring_number).KD=kd; % damping coefficient

    particle = particle +1;
  
end

%%

for i=1:length(spring)
    LINJE_a1(i)=line([x(spring(i).from,1),x(spring(i).to,1)],[x(spring(i).from,2),x(spring(i).to,2)]);
end

Angularmomentum = zeros(1,N);
Springlength = zeros(1,N);
Angularfreq = zeros(1,N);
Es = zeros(1,N);
Ek = zeros(1,N);
Ep = zeros(1,N);

F = zeros(NP,dim);

%%

for k = 1:NS

    r = x(spring(k).from,:)-x(spring(k).to,:);
    rdot = v(spring(k).from,:)-v(spring(k).to,:);

    f = -((spring(k).KS)*(norm(r)-spring(k).length)+spring(k).KD*(dot(rdot,r)/norm(r))).*r/norm(r);

    F(spring(k).from,:) = F(spring(k).from,:) + f;
    F(spring(k).to,:) = F(spring(k).to,:) - f;

    Es = Es + spring(k).KS*(norm(r)-spring(k).length).^2./2;
end

v = v - dt*inv(M)*F;%./2;

%%

for n = 1:N
    
    F = zeros(NP,dim);
    es = 0;

    for k = 1:NS

        r = x(spring(k).from,:)-x(spring(k).to,:);
        rdot = v(spring(k).from,:)-v(spring(k).to,:);

        f = -((spring(k).KS)*(norm(r)-spring(k).length)+spring(k).KD*(dot(rdot,r)/norm(r))).*r/norm(r);

        F(spring(k).from,:) = F(spring(k).from,:) + f;
        F(spring(k).to,:) = F(spring(k).to,:) - f;

        es = es + spring(k).KS*(norm(r)-spring(k).length).^2./2;
    end
    F(:,dim) = F(:,dim)-m'*g;
    v_prev = v;
    v = v + dt*inv(M)*F;
    x = x + dt*v;
    
    springlength(n) = norm(r);
    Angularmomentum(n) = m*(x(:,1).*v(:,2)-x(:,2).*v(:,1));
    Angularfreq(n) = norm(v(1,:)./(norm(r)));
    Ek(n) = sum(m.*(sum(((v+v_prev)./2).^2)./2));
    Ep(n) = sum(g.*m.*x(:,2)');

    %uppdaterad grafik
    for i = 1:NS
        set(LINJE_a1(i),'XData',[x(spring(i).from,1),x(spring(i).to,1)],'YData',[x(spring(i).from,2),x(spring(i).to,2)]);
    end
    for i=1:NP
        set(BOLL_a(i),'Position',[x(i,1)-Ra,x(i,2)-Ra,2*Ra,2*Ra]);
    end
    Es(n) = es;
    pause(dt/1000);
end

%%
Etot = Es + Ek + Ep;


for i = 1:1:length(time)
    energy = Etot(i)/Etot(1);
    if energy <= 0.1
        break;
    end
end
TimeTo10percent = time(i)


figure(2)
plot(time,Etot)
hold on
plot(time,Ek)
hold on
plot(time,Ep)
hold on
plot(time,Es)
hold off

figure(3)
plot(time,Angularmomentum)
hold off

figure(4)
plot(time,Angularfreq)
hold on
plot(time,springlength)