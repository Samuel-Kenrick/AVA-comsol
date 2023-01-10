% particles and positions
NP = 3;
M = [1,1,1]
L = [1,1,1]
x0 = [-L/2,0;L/2,0;0,L/2]

% velocity
vini = 0.01;
v0 = [vini,0;0,-vini;-vini/2,-vini/2]


%timestep
timesteps = 1500
dt = 0.1
time = 1:dt:timesteps*dt
N = length(time)

m = 1



for n = 4:N

    F(n) = 0.5*(0.5-rand(NP,2))-0.1*v(n)
    v(n+1/2)=v(n-1/2)+dt*inv(M)*F(n)
    X(n+1)=x(n)+dt*v(n+1/2)
end
    
