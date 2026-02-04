%Code for ACS Perspective : Considerations for Gas- and Liquid-Fed Electrocatalytic Coupling: A Case Study with Carbon Dioxide and Nitrate
%Sahil Halarnkar^1, Andrew Ruba^2, Nirala Singh*1,2
%* Corresponding Author: snirala@umich.edu
%^ These authors contributed equally.
%1. Department of Chemical Engineering, University of Michigan, Ann Arbor, MI, 
%2. Macromolecular Science & Engineering, University of Michigan, Ann Arbor, MI
%% Kinetics Diffusion model
%% Parameters 
k1=1;
k2=1;
k3=1;
k4=1;
L=1;
cA=1;
cB=1;
cC=1;
Da=1;
Db=1;
Dc=1;
TA1=sqrt(k1*L*L/Da);
TA3=sqrt(k3*L*L/Da);
TB2=sqrt(k2*L*L/Db);
TB3=sqrt(k3*L*L/Db);
TC1=sqrt(k1*L*L/Dc);
TC2=sqrt(k2*L*L/Dc);
TC3=sqrt(k3*L*L/Dc);
TC4=sqrt(k4*L*L/Dc);
% Parameters
a1=TA1^2*cC;
a2=TA3^2*cC*cB;
a3=TB2^2*cC;
a4=TB3^2*cC*cA;
a5=TC1^2*cA;
a6=TC2^2*cB;
a7=TC3^2*cA*cB;
%% Figure 2b
% Initial mesh and guess for solution
xmesh = linspace(0, 1, 20);
solinit = bvpinit(xmesh, [1 0 1 0 1 0]);

% Solve BVP
sol = bvp4c(@(x,Y) GDE(x,Y,a1,a2,a3,a4,a5,a6,a7), @bcfun, solinit);

% Extract solution
x = sol.x;
y = sol.y(1,:)*cA;
z = sol.y(3,:)*cB;
v = sol.y(5,:)*cC;

figure;
plot(x,y,'r-','LineWidth',5); hold on;
plot(x,z,'b--','LineWidth',5); hold on;
plot(x,v,'k-.');
legend('N','C','H');
xlabel('x');
ylim([0 1]);
%% Figure S1b

% Initial mesh and guess for solution
xmesh = linspace(0, 1, 20);
solinit = bvpinit(xmesh, [1 0 1 0 1 0]);

% Solve BVP
sol = bvp4c(@(x,Y) GDE(x,Y,a1,a2,a3,a4,a5,a6,a7), @bcfun_opp, solinit);

% Extract solution
x = sol.x;
y = sol.y(1,:)*cA;
z = sol.y(3,:)*cB;
v=sol.y(5,:)*cC;
figure;
plot(x,y,'r-','LineWidth',5); hold on;
plot(x,z,'b--','LineWidth',5); hold on;
plot(x,v,'k-.');
legend('N','C','H');
xlabel('x');
ylim([0 1]);


%% Figure 2 c-e

Li=logspace(-2,1,100); %Choose range 
i=1;
eta=0; % Internal effectiveness factor
select_CNW=0; % Coupling reaction selectivity
rr_CNW=0; % Coupling reaction rate 

% Choose parameter to be varied as loop variable 
for L=Li
    
    % Initial mesh and guess for solution
    xmesh = linspace(0, 1, 20);
    solinit = bvpinit(xmesh, [1 0 1 0 1 0]);
    
    % Solve BVP
    sol = bvp4c(@(x,Y) GDE(x,Y,a1,a2,a3,a4,a5,a6,a7), @bcfun, solinit);
    
    % Extract solution
    x = sol.x;
    y = sol.y(1,:)*cA;
    z = sol.y(3,:)*cB;
    v = sol.y(5,:)*cC;
    
    eta(i)=(trapz(x,y.*z.*v))/(cA*cB*cC);
    select_CNW(i)=(k3*trapz(x,y.*v.*z)/(k1*trapz(x,y.*v)+k2*trapz(x,z.*v)+k3*trapz(x,y.*v.*z)));
    rr_CNW(i)=k3*trapz(x,y.*v.*z);
    i=i+1;
end
figure
plot(Li,eta);
hold on
%% Selectivity conversion model

%% Parameters

% Rate constants
k.kN1  = 0.1;
k.kN2  = 0.1;
k.kC1  = 0.1;
k.kC2  = 0.1;
k.kCN1 = 0.1;
k.kCN2 = 0;
k.kCN3 = 0;
k.kCN4 = 0;

% Initial concentrations
C0 = [
    1.0;  % NO3-
    1.0;  % CO2
    0.0;  % NO2-
    0.0;  % CO
    0.0   % CN
];

% Time span
tspan = [0 20];
%% Figures S3-S6
% Solve ODEs
[t, C] = ode45(@(t,C) kinetics_CN(t, C, k), tspan, C0);

% Plot
figure;
X=1-C(:,1); % Conversion of NO3-
plot(t,C(:,1),'r','LineWidth',5)
hold on;
plot(t,C(:,2),'b--','LineWidth',5)
hold on;
plot(t,C(:,3),'Color',[1 0.5 0],'LineWidth',5)
hold on;
plot(t,C(:,4),'LineStyle','--','Color',[0.75 0.25 0.5],'LineWidth',5)
hold on;
plot(t,C(:,5),'LineStyle',":",'Color','k','LineWidth',5)
hold on;
legend('NO_3^-','CO_2','NO_2^-','CO','CN')
xlim([0 15])
xlabel('Space Time');
ylabel('Concentration');
%% Figure 3b

% Time span
tspan = [0 20];
k.kCN1 = 0.1;
k.kCN2 = 0;
k.kCN3 = 0;
k.kCN4 = 0;
[t, C1] = ode45(@(t,C) kinetics_CN(t, C, k), tspan, C0);
k.CN1=0;
k.kCN2 = 0.1;
k.kCN3 = 0;
k.kCN4 = 0;
[t, C2] = ode45(@(t,C) kinetics_CN(t, C, k), tspan, C0);

k.kCN1 = 0;
k.kCN2 = 0;
k.kCN3 = 0.1;
k.kCN4 = 0;
[t, C3] = ode45(@(t,C) kinetics_CN(t, C, k), tspan, C0);

k.kCN1 = 0;
k.kCN2 = 0;
k.kCN3 = 0;
k.kCN4 = 0.1;
[t, C4] = ode45(@(t,C) kinetics_CN(t, C, k), tspan, C0);
X1=1-C1(:,1);
X2=1-C2(:,1);
X3=1-C3(:,1);
X4=1-C4(:,1);
plot(X1,(C1(:,1).*C1(:,2))./(C1(:,1)+C1(:,2)))
hold on;
plot(X2,C2(:,1).*C2(:,4)./(C2(:,1)+C2(:,2)))
hold on
plot(X3,C3(:,3).*C3(:,2)./(C3(:,1)+C3(:,2)))
hold on
plot(X4,C4(:,3).*C4(:,4)./(C4(:,1)+C4(:,2)))
hold off
legend('1','2','3','4')
%% Functions
% ---- System of ODEs ----
function dYdx = GDE(x,Y,a1,a2,a3,a4,a5,a6,a7)
    y  = Y(1); yp = Y(2);
    z  = Y(3); zp = Y(4);
    v=Y(5); vp=Y(6);

    dYdx = zeros(6,1);
    dYdx(1) = yp;               % y'
    dYdx(2) = a1*y*v + a2*y*z*v;    % y''
    dYdx(3) = zp;               % z'
    dYdx(4) = a3*z*v + a4*y*z*v;
    dYdx(5) = vp;               % z'
    dYdx(6) = a5*z*v + a6*y*v + a7*y*z*v; % z'' %remove a8
end

% ---- Boundary conditions ----
function res = bcfun(Y0,YL)
    % Boundary conditions:
    % y(0)=1, y'(L)=0, z(L)=1, z'(0)=0
    res = [Y0(1) - 1;    % y(0) = 1
           YL(2) - 0;    % y'(L) = 0
           YL(3) - 1;    % z(L) = 1
           Y0(4) - 0     % z'(0) = 0
           Y0(5)-1;
           YL(6)-0;];   
end
function res = bcfun_opp(Y0,YL)
    % Boundary conditions:
    % y(0)=1, y'(L)=0, z(L)=1, z'(0)=0
    res = [Y0(1) - 1;    % y(0) = 1
           YL(2) - 0;    % y'(L) = 0
           Y0(3) - 1;    % z(L) = 1
           YL(4) - 0     % z'(0) = 0
           Y0(5)-1;
           YL(6)-0;];   
end

% ---- Initial guess ----
function g = guess(x)
    % Guess for [y, y', z, z']
    g = [1 - x/10;   % y ~ linear decay
         -0.1;       % y' ~ small slope
          x/10;      % z ~ linear growth
          0.1];      % z' ~ small slope
end
function dCdt = kinetics_CN(t, C, k)
% C = [NO3; CO2; NO2; CO; CN]
% k = struct of rate constants
dCdt=zeros(5,1);
NO3 = C(1);
CO2 = C(2);
NO2 = C(3);
CO  = C(4);
CN  = C(5); 

dCdt(1) = -k.kN1*NO3 -k.kCN1*NO3*CO2 -k.kCN2*NO3*CO;

dCdt(2) = -k.kC1*CO2 -k.kCN1*NO3*CO2 -k.kCN3*NO2*CO2;

dCdt(3) =  k.kN1*NO3  -k.kN2*NO2  -k.kCN3*NO2*CO2  -k.kCN4*NO2*CO;

dCdt(4) =  k.kC1*CO2  -k.kC2*CO  -k.kCN2*NO3*CO  -k.kCN4*NO2*CO;

dCdt(5) =  k.kCN1*NO3*CO2  +k.kCN2*NO3*CO  +k.kCN3*NO2*CO2  +k.kCN4*NO2*CO;


end
