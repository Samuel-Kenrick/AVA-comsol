% positions
L = 1
x(1) = [-L/2,0]
x(2) = [L/2,0]
x(3) = [0,L/2]

% velocity
vini = 0.01
v(1) = [vini,0]
v(2) = [0,-vini]
v(3) = [-vini/2,-vini/2]

%timestep
timesteps = 1500
dt = 0.1
time = 1:0.1:timesteps
N = length(time)

m = 1



for i = 2:N
    
    v()
F = 0.5*(0.5-rand(NP,2))-0.1*V
