% TGA simulation

%% variable initalization

tic;

clear all vars 
global ycoeff afac nfac ea istart g_index s_index  MW rhos0 gsp nsp  masslossrate 

load ('ranzi_pyro_kinetics_gentile2017.mat');
MW(47) = 28;
ycoeff(47,:) = 0;
g_index = [3 4 5 6 7 8 9 10 11 12 13 14 16 20 21 22 29 30 31 33 34 35 47];
gsp = length(g_index);
nsp = 47;
s_index = [1 2 15 17 18 19 23 24 25 26 27 28 32 36 37 38 39 40 41 42 43 44 45 46];

%% sizing and array initialization

Mesh.Jnodes = 1;
sample_height = 5e-3;
Mesh.z = linspace(0,sample_height,Mesh.Jnodes);
Mesh.dz = sample_height/(Mesh.Jnodes);
Mesh.a = 5e-3; 
Mesh.dv = Mesh.a * Mesh.dz;

yprime0 = zeros(Mesh.Jnodes*nsp+2*Mesh.Jnodes,1);
rhos0 = zeros(Mesh.Jnodes,1);
rhos_mass0 = zeros(Mesh.Jnodes,1);

%% set initial composition

m0 = zeros(47,Mesh.Jnodes);

m0(1,:) = (0.24092)/MW(1); % CELL 
m0(17,:) = (0.156317)/MW(17); % HCE
m0(24,:) = (0.13672/3)/MW(24); % LIGH
m0(25,:) = (0.13672/3)/MW(25); % LIGO
m0(23,:) = (0.13672/3)/MW(23); % LIGC
m0(38,:) = (0.023269)/MW(38); % TGL
m0(37,:) = (0.003384)/MW(37); % TANN
m0(39,:) = (0.30)/MW(39); % moisture

%% initialize mass calcs 

mass0 = m0.*MW; % kg
% yi0 = mass0(s_index,1)./100;
sample_density = 492; % [kg/m3] 400
rhos_mass0 = rhos_mass0+sample_density; % +100
sample_mass = Mesh.a*sample_height*rhos_mass0(1); % 48.2817e-6; 
mass0 = mass0./sum(mass0(:,1))*sample_mass./Mesh.Jnodes;
y0 = [mass0(:); rhos_mass0(:)];

%% initial conditions

T0 = 318; % initial temperature
Tend = 1200; % final temperature
dt = 1;
beta = 50/60;  % rate of temperature change [K/s]
nstep = fix((Tend-T0)/beta)*(1/dt);
time = 0;
t = zeros(nstep+1,1); 
yy = zeros(nstep+1,length(y0)); 
t(1)= 0;
yy(1,:) = y0;
ye = zeros(length(t),length(g_index));
j0 = zeros(length(t),1);
mlr = zeros(length(t),1);
T = zeros(length(t),1); T(1) = 318; % adjust from 300
Mg = MW(g_index)*1e-3;

options = odeset('RelTol',1.e-4,'AbsTol',1e-5, 'NonNegative', 1, 'BDF',0, 'MaxOrder',2);

%% time integration 

for i=1:nstep
    tspan = [t(i) t(i)+dt];
    [t2,a] = ode113(@(t,y)yprime(time,y,Mesh,T(i)),tspan,yy(i,:),options);
    temp = a(end,:);
    temp(temp<0)=1e-30;
    mlr(i+1) = masslossrate;
    yy(i+1,:) = temp;
    T(i+1) = T(i)+ beta*dt;
    t(i+1) = t(i) + dt;    
end

%% plot

normalized_mass = yy(:,end)/sample_density; 

figure(1); clf
hold on;
plot(T, normalized_mass); %yy(:,end)/sample_density); % normalize by sample density
% xlim([300 1250]);
% ylim([0 100]);
xlabel('Temp [K]');
ylabel('mass %');
title('mass % evolution wrt T');

figure(2); clf
plot(T, -mlr);
xlabel('Temperature [K]');
ylabel('mlr, DTG');
title('Mass loss rate (mlr, DTG) wrt T');
hold off;

toc;
runtimetga = toc;

%% define functions 

function [dydt] = yprime(t,yy,Mesh,T)

global ycoeff afac nfac ea istart s_index g_index MW nsp masslossrate yje

    wdot_mass = zeros(nsp,Mesh.Jnodes);
    k = zeros(28,Mesh.Jnodes);
    m = zeros(nsp,Mesh.Jnodes);
    rho_s_mass = zeros(Mesh.Jnodes,1);
    drhosdt = zeros(Mesh.Jnodes,1);
    mprime = zeros(nsp,Mesh.Jnodes);
    yje = zeros(length(g_index),1);
    
    for i=1:Mesh.Jnodes
        temp=yy(nsp*(i-1)+1:nsp*(i-1)+nsp);
        temp(temp<0)=1e-30;
        m(:,i)=temp;
        m(:,i)= m(:,i)./MW;
    end
    
    yi = zeros(length(s_index),Mesh.Jnodes);
    rho_s_mass(:) = yy(nsp*Mesh.Jnodes+1:end);
    
    R = 8.314; 
    
    for i=1:Mesh.Jnodes
        yi(:,i) = m(s_index,i).*MW(s_index)./sum(m(s_index,i).*MW(s_index));
        k(:,i) = afac .*((T(i)).^nfac).* exp(-ea ./(R*T(i)));
        mprime(:,i) = ycoeff*(k(:,i).*m(istart,i)).*MW; % dm/dt
        wdot_mass(:,i) = mprime(:,i)./ Mesh.dv;
    end 
  
    for i=1:Mesh.Jnodes
        drhosdt(i) = - sum(wdot_mass(g_index,i));
    end

    masslossrate = sum(drhosdt);  
    dydt = [mprime(:); drhosdt(:)];  
end

function phi = phii(yi,rho_s_mass)
    
    global MW s_index
    
    s_density = [9.37000000000000;9.37000000000000;25;9.87000000000000;11.5000000000000;11.5000000000000;...
        11.5000000000000;12.1200000000000;5.88000000000000;3.48000000000000;3.59000000000000;...
        5.88000000000000;4;7.29000000000000;5.76000000000000;7.22000000000000;5;1.67000000000000;...
        55;0.00369448575008421;0.00580475748165167;0.00541504238833635;0.0806563778419628;...
        0.0101350128633761;0.00507436386823035;0.00579578562602859]; 

    phi = 1-sum(yi./(s_density.*MW(s_index)))*rho_s_mass;
end

function rho_sm = rhos_mass(yi,phi)
% kg/m3
    global MW s_index
    
    s_density = [9.3745;9.3745;25;9.87000000000000;11.5050;11.5050;...
        11.5050;12.1200000000000;5.8852;3.4826;3.59000000000000;...
        5.88000000000000;4;7.29000000000000;5.76000000000000;7.22000000000000;5;1.67000000000000;...
        55;0.00369448575008421;0.00580475748165167;0.00541504238833635;0.0806563778419628;...
        0.0101350128633761;0.00507436386823035;0.00579578562602859];

    rho_sm = (1-phi)/sum(yi./(s_density.*MW(s_index)));
end 