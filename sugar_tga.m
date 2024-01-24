% tga sim for sugars only
% by DIBA

%% initialize 
clear all vars
global ycoeff afac nfac ea istart g_index s_index MW gsp nsp_len masslossrate 
species = {'sugar','sugar1','sugar2','taro1','taro2','h2o','taro3','char','gco2','gch4','gc2h4','gco','gcoh2'};

ycoeff=[-1	0	0
        0.47 -1	0
        0.53 0	-1
        0	0.6	0
        0	0.4	0
        0	0.4	0.88
        0	0	0.13
        0	0	1.6
        0	0.9	1.3
        0	0	0.25
        0	0	0.1
        0	0	0.62
        0	0	0.73];

% check if TARO species end up in gas phase
g_index=[4 5 6 7]; % 4?5? try yes then no
s_index=[1 2 3 8 9 10 11 12 13];
gsp = length(g_index);
nsp = [g_index, s_index]; % all species

istart=[1 2 3];

% add molecular weights
MW = [176; 176; 176; 176; 114; 18; 290; 12; 44; 16; 28; 28; 30];

afac = [.8e10 .15e5 .2e2];
nfac = [0 0 0];

% activation energy, convert to J/mol
ea = [26000 16000 20000]*4.184;

Mesh.Jnodes = 1; % should be 1 for tga
sample_height = 1e-2;  % m?
Mesh.z = linspace(0,sample_height,Mesh.Jnodes);
Mesh.dz = sample_height/(Mesh.Jnodes);
Mesh.a = 1e-2^2;
Mesh.dv = Mesh.a * Mesh.dz;

nsp_len = length(nsp); % total (gas and solid) species array length
yprime0 = zeros(Mesh.Jnodes*nsp_len+2*Mesh.Jnodes ,1);
rhos_mass0 = zeros(Mesh.Jnodes,1);

% set initial composition
m0 = zeros(13,Mesh.Jnodes);
m0(1) = 1;  % line added to specify only reactant at start time (first 
% reaction - sugar)

% define # moles ea. species in question (sugar reactions)

rhos_mass0 = rhos_mass0+100;

sample_mass = Mesh.a*sample_height*rhos_mass0(1);
m0 = m0./sum(m0(s_index,1))*sample_mass./Mesh.Jnodes; %%%%%%%%%%%%%%%%%%

y0 = [m0(:); rhos_mass0(:)];

% specify initial conditions
T0 = 300; % initial temperature
Tend = 700; % final temperature. is there a way to not pre-set this?
dt = 1;
beta = 10/60; %rate of temperature change (K/s)
nstep = fix((Tend-T0)/beta)*5; %*100
time = 0;
t = zeros(nstep+1,1); 
yy = zeros(nstep+1,length(y0)); 
t(1)= 0;
yy(1,:) = y0;
ye = zeros(length(t),length(g_index));
j0 = zeros(length(t),1);
mlr = zeros(length(t),1);
T = zeros(length(t),1); T(1) = 300;

options = odeset('RelTol',1.e-4,'AbsTol',1e-5, 'NonNegative', 1, 'BDF',0, 'MaxOrder',2);

%% begin iterating through mesh
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
figure(1); clf
hold on;
plot(T, yy(:,end));
xlabel('Temp [K]');
ylabel('mass loss wrt T');
title('Mass lost wrt t (mlr)');

figure(2); clf
plot(T, mlr);
xlabel('Temperature [K]');
ylabel('mlr, DTG');
title('Mass loss rate (mlr, DTG) wrt T');
hold off;

%% define functions
function [dydt] = yprime(t,yy,Mesh,T)

global ycoeff afac nfac ea istart s_index g_index MW nsp_len masslossrate yje

    wdot_mass = zeros(nsp_len,Mesh.Jnodes);
    k = zeros(3,Mesh.Jnodes);
    m = zeros(nsp_len,Mesh.Jnodes);
    drhosdt = zeros(Mesh.Jnodes,1);
    mprime = zeros(nsp_len,Mesh.Jnodes);
    
    yje = zeros(length(g_index),1);
    
    for i=1:Mesh.Jnodes
        temp=yy(nsp_len*(i-1)+1:nsp_len*(i-1)+nsp_len);
        temp(temp<0)=1e-30;
        m(:,i)=temp;
        %m(:,i)= transpose(m)./MW;
        m(:,i)= m(:,i)./MW;
    end
    
    yi = zeros(length(s_index),Mesh.Jnodes);
   
    
    R = 8.314; 
    
    for i=1:Mesh.Jnodes
    
        yi(:,i) = m(s_index,i).*MW(s_index)./sum(m(s_index,i).*MW(s_index));
        k(:,i) = afac .*((T(i)).^nfac).* exp(-ea ./(R*T(i)));
        mprime(:,i) = ycoeff*(k(:,i).*m(istart,i)).*MW; %dmdt
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

    global MW s_index
    
    s_density = [9.3745;9.3745;25;9.87000000000000;11.5050;11.5050;...
        11.5050;12.1200000000000;5.8852;3.4826;3.59000000000000;...
        5.88000000000000;4;7.29000000000000;5.76000000000000;7.22000000000000;5;1.67000000000000;...
        55;0.00369448575008421;0.00580475748165167;0.00541504238833635;0.0806563778419628;...
        0.0101350128633761;0.00507436386823035;0.00579578562602859];

    rho_sm = (1-phi)/sum(yi./(s_density.*MW(s_index)));
end