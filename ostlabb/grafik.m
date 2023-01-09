function grafik()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exempel av spinnande boll, boll som deformeras och
% roterande rektangel som ror sig.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametervaerden for spinnande bollen
xa=2;    % Start-position (centrum)
ya=5;
Da=0;    % Vinkeln
Ra=0.7;  % Radie
La=0.4;  % Kryssets storlek

% Parametervaerden for deformerad boll
xb=4;    % Start-position (centrum)
yb=3;
Rb=0.7;  % Radie
du=0;    % Ovre intryckning
dd=0;    % Undre intryckning

% Parametervaerden for rektangel
xc=6;    % Start-position (centrum)
yc=2;
Dc=0;    % Vinkeln som rektangeln aer vriden (moturs)
Bc=1;    % Basen
Hc=2;    % Hojden

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utritning av spinnande bollen och dess kryss
% Bollens center = [xa ,ya]
% Rita en rektangel med start i bollens center minus bollens radie i x och
% y, med bredd och hojd av dubbla radien.
% Sedan saetts kurvaturen till 1 (100%) i respektive riktning och gor att
% kvadraten (rektangeln) blir en cirkel.
% Bollen:
BOLL_a=rectangle('Position',[xa-Ra,ya-Ra,2*Ra,2*Ra],...
    'Curvature',[1,1],'EdgeColor','r','FaceColor','w');
% Krysset (tvaa korsade linjer):
xk1_1=-La/2*sin(Da);    % Aendpunkter till linjerna som
yk1_1=La/2*cos(Da);     % utgor kryss nummer 1
xk1_2=La/2*sin(Da);
yk1_2=-La/2*cos(Da);
xk2_1= La/2*cos(Da);    % Aendpunkter till linjerna som
yk2_1=La/2*sin(Da);     % utgor kryss nummer 2
xk2_2=-La/2*cos(Da);
yk2_2=-La/2*sin(Da);
% Rita linjerna
LINJE_a1=line([xk1_1+xa,xk1_2+xa],[yk1_1+ya,yk1_2+ya]);
LINJE_a2=line([xk2_1+xa,xk2_2+xa],[yk2_1+ya,yk2_2+ya]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utritning av deformerad boll
BOLL_b=rectangle('Position',...
    [xb-(Rb+du+dd),yb-(Rb-dd),2*(Rb+du+dd),2*(Rb-du-dd)],...
        'Curvature',[1,1],'EdgeColor','b','FaceColor','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utritning av rektangel
% Beraekna forst hornens positioner daa den
% inte aer vriden, sett fraan rektangelns centrum
x1c=Bc/2;    % Ovre hoger horn (1)
y1c=Hc/2;
x2c=-Bc/2;   % Ovre vaenster horn (2)
y2c=Hc/2;
x3c=-Bc/2;   % Nedre vaenster horn (3)
y3c=-Hc/2;
x4c=Bc/2;    % Nedre hoger horn (4)
y4c=-Hc/2;
% Vrid dessa punkter en vinkel Dc (moturs)
X1c=cos(Dc)*x1c-sin(Dc)*y1c;   % Horn 1
Y1c=sin(Dc)*x1c+cos(Dc)*y1c;
X2c=cos(Dc)*x2c-sin(Dc)*y2c;   % Horn 2
Y2c=sin(Dc)*x2c+cos(Dc)*y2c;
X3c=cos(Dc)*x3c-sin(Dc)*y3c;   % Horn 3
Y3c=sin(Dc)*x3c+cos(Dc)*y3c;
X4c=cos(Dc)*x4c-sin(Dc)*y4c;   % Horn 4
Y4c=sin(Dc)*x4c+cos(Dc)*y4c;
% Rita rektangeln som segment av linjer
REKTANGEL_c=line([X1c+xc,X2c+xc,X3c+xc,X4c+xc,X1c+xc],...
    [Y1c+yc,Y2c+yc,Y3c+yc,Y4c+yc,Y1c+yc]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instaellningar for grafik-fonstret
daspect([1,1,1]);   % Skalar axlarna lika
axis([0,8,0,6]);    % Axlarnas intervall
hold on             % Laas fast dessa instaellningar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Huvudslingan
for n=1:50

    % Updatera positioner, vinkel, intryckning
    ya=ya-0.05;
    Da=Da+0.075;
    du=0.2*Rb*sin(0.3*n);
    dd=0.2*Rb*sin(0.3*n);
    yc=yc+0.05;
    Dc=Dc+0.15;

    % Uppdatera grafiken:

    % Spinnande bollen och krysset
    % Uppdatera bollens position med ''set''
    set(BOLL_a,'Position',[xa-Ra,ya-Ra,2*Ra,2*Ra]);
    % De nya linjernas aendpunkter
    xk1_1=-La/2*sin(Da);   % Kryss nummer 1
    yk1_1=La/2*cos(Da);
    xk1_2=La/2*sin(Da);
    yk1_2=-La/2*cos(Da);
    xk2_1= La/2*cos(Da);   % Kryss nummer 2
    yk2_1=La/2*sin(Da);
    xk2_2=-La/2*cos(Da);
    yk2_2=-La/2*sin(Da);
    % Uppdatera linjerna med ''set''
    set(LINJE_a1,'XData',[xk1_1+xa,xk1_2+xa],...
        'YData',[yk1_1+ya,yk1_2+ya]);
    set(LINJE_a2,'XData',[xk2_1+xa,xk2_2+xa],...
        'YData',[yk2_1+ya,yk2_2+ya]);

    % Deformerade bollen
    % Uppdatera bollens position och
    % intryckning med ''set''
    set(BOLL_b,'Position',[xb-(Rb+du+dd),yb-(Rb-dd),...
        2*(Rb+du+dd),2*Rb-du-dd]);

    % Rektangeln
    % Beraekna hornens nya positioner
    X1c=cos(Dc)*x1c-sin(Dc)*y1c;   % Horn 1
    Y1c=sin(Dc)*x1c+cos(Dc)*y1c;
    X2c=cos(Dc)*x2c-sin(Dc)*y2c;   % Horn 2
    Y2c=sin(Dc)*x2c+cos(Dc)*y2c;
    X3c=cos(Dc)*x3c-sin(Dc)*y3c;   % Horn 3
    Y3c=sin(Dc)*x3c+cos(Dc)*y3c;
    X4c=cos(Dc)*x4c-sin(Dc)*y4c;   % Horn 4
    Y4c=sin(Dc)*x4c+cos(Dc)*y4c;
    % Uppdatera rektangeln med ''set''
    set(REKTANGEL_c,'XData',[X1c+xc,X2c+xc,X3c+xc,X4c+xc,X1c+xc],...
        'YData',[Y1c+yc,Y2c+yc,Y3c+yc,Y4c+yc,Y1c+yc]);

    pause(0.1);    % En paus for dra ner paa tempot
end
