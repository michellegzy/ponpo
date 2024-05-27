%% 1D pyrolysis solver by DBehnoudfar and Michelle Gee

tic
clear vars

global reactions afac nfac ea istart qs g_index s_index MW gsp ssp tempflux p0 Kd yj0
%% load species data and kinetics parameters

load ('updated_kinetics3.mat');

% reactions(49,:) = 0; % reaction coefficients for N2
% reactions(50,:) = 0; % reaction coefficients for ASH

% gas-phase species indices
%g_index = [3 4 5 6 7 8 9 10 11 12 13 14 16 20 21 22 29 30 31 33 34 35 47 48];
gsp = length(g_index); % # of gas-phase species

% solid-phase species indices
%s_index = [1 2 15 17 18 19 23 24 25 26 27 28 32 36 37 38 39 40 41 42 43 44 45 46 50];
ssp = length(s_index); % # of solid-phase species
MW = MW * 1e-3; % conversion from g/mol to kg/mol
solid_densities = s_density.*MW(s_index)*1000; % convert array from mol/cm3 to g/cm3

%% setup mesh 

Mesh.Jnodes = 3; % mesh size
sample_height = 3.8e-2; % [m]
Mesh.dz = sample_height/(Mesh.Jnodes);
Mesh.a = 1e-2^2; % surface area, [m2]
Mesh.dv = Mesh.a * Mesh.dz;

%% initialize variables

T0 = zeros(Mesh.Jnodes,1);
yprime0 = zeros(Mesh.Jnodes*ssp+2*Mesh.Jnodes,1);
rhos_mass0 = zeros(Mesh.Jnodes,1); % density
Tinitial = 300; % [K]

for j = 1:Mesh.Jnodes
    T0(j) = Tinitial; 
end

m0 = zeros((ssp+gsp),Mesh.Jnodes); % mole storage matrix

%% set initial composition and conditions

% define # moles at each node (column) for corresp. species (row)
m0(1,:) = 0.4254/MW(1); % CELL 
m0(17,:) = 0.1927/MW(17); % HCE
m0(24,:) = 0.0998/MW(24); % LIGH
m0(25,:) = 0.0482/MW(25); % LIGO
m0(23,:) = 0.1658/MW(23); % LIGC
m0(38,:) = 0.0326/MW(38); % TGL
m0(37,:) = 0.0354/MW(37); % CTANN
m0(39,:) = 0.05/MW(39); % moisture
m0(50,:) = 0.05/MW(49); % ash

mass0 = m0.*MW; % [kg]
yi0 = mass0(s_index,1)./sum(mass0(s_index,1));
sample_density = 380; % [kg/m3]
rhos_mass0 = rhos_mass0+sample_density;
sample_mass = Mesh.a*sample_height*rhos_mass0(1); 
mass0 = mass0./sum(mass0(s_index,1))*sample_mass./Mesh.Jnodes;

% set environment
p0 = 1.013e5; % pressure [Pa]
yj0 = zeros(gsp,1); 
R = 8.314; % univ gas ct [J/molK]

% the next 2 lines define the gas environment for computations by the last
% 2 species in the gas index array
yj0(end-1) = 0.21; % O2, gas-phase species mass fraction
yj0(end) = 0.79; % N2
M = 1/sum(yj0./MW(g_index)); 

rhog0 = zeros(Mesh.Jnodes,1); % gas-phase density
rhog0 = rhog0 + (p0)*M/(R*T0(1));
rhogphi0 = rhog0*phii(yi0,rhos_mass0(1), solid_densities); % added input arg solid_densities
rgpy0 = zeros(gsp,Mesh.Jnodes); % rhog * phi * y
rgpy00 = rhogphi0(1).*yj0; 
for i = 1:Mesh.Jnodes
    rgpy0(:,i)= rgpy00;
end
rgpy0 = reshape(rgpy0,gsp*Mesh.Jnodes,1);

qs = 40000; % input heat flux [w/m2]
y0 = [rhogphi0(:); rgpy0(:)]; 
y10 = [mass0(:); T0(:); rhos_mass0(:)];
phi = phii(yi0,rhos_mass0(1), solid_densities); % fuel porosity % added input arg solid_densities
Kd = 1e-10; % porous fuel permeability

%% ode solver options

dt = .1; 
nstep = 10; 
time = 0; 
t = zeros(nstep+1,1); 
yy = zeros(nstep+1,length(y0)); % species transport equation solution matrix
yy1 = zeros(nstep+1,length(y10)); % heat equation solution matrix

t(1) = 0;
yy(1,:) = y0;
yy1(1,:) = y10;
ye = zeros(length(t),length(g_index));
j0 = zeros(length(t),1);
Ts = zeros(length(t),1); 
Ts(1) = 300;

options1 = odeset('RelTol',1e-4,'AbsTol',1e-5,'BDF',0,'MaxOrder',1);

%% time integration 

for i = 1:nstep
    tspan = [t(i) t(i)+dt];
    [~,b] = ode15s(@(t,y)yprime1(time,y,Mesh),tspan,yy1(i,:),options1);
    [~,a] = ode15s(@(t,y)yprime(time,y,Mesh,yy1(i,:)),tspan,yy(i,:),options1);

    temp = a(end,:);
    temp(temp<0) = 1e-30;
    j0(i+1) = tempflux;
	tempflux; % comment out? 
    yje = temp(:,end-gsp+1:end)./sum(temp(:,end-gsp+1:end),2);    
    ye(i+1,:) = yje;
    yy(i+1,:) = temp;
    yy1(i+1,:) = b(end,:);
    Ts(i+1) = yy1(i+1,ssp*Mesh.Jnodes+Mesh.Jnodes);
    t(i+1) = t(i) + dt;
end

%% plotting and related operations

dimensionless_rho = yy(:,end)/yy(1,end);

% figure(1); clf
% hold on;
% plot(t, dimensionless_rho);
% xlabel('time [s]');
% ylabel('density ratio');
% title('initial / final density over time');

% figure(2); clf
% plot(T, yy(:,end));
% %xlim([300 1250]);
% %ylim([0 100]);
% xlabel('Temp [K]');
% ylabel('mass %');
% title('mass % evolution wrt T');
 
% figure(3); clf
% plot(T, -mlr);
% xlabel('Temperature [K]');
% ylabel('mlr, DTG');
% title('Mass loss rate (mlr, DTG) wrt T');
% hold off;

toc
%% define functions 

% ODE function for species trassport equation
function [dydt] = yprime(t,yy,Mesh,yy1, reactions, afac, nfac, ea, istart, s_index, g_index, MW, gsp, ssp, tempflux, p0, Kd)

    %global reactions afac nfac ea istart s_index g_index MW gsp ssp tempflux p0 Kd

    wdot_mass = zeros(ssp,Mesh.Jnodes); % species production rate
    %k = zeros(28,Mesh.Jnodes); % reaction rate coefficient
    k = zeros(size(reactions,2), Mesh.Jnodes); % reaction rate coefficient
    m = zeros(ssp,Mesh.Jnodes); % mass
    phi = zeros(Mesh.Jnodes,1); % porosity
    kb = zeros(Mesh.Jnodes,1); % conductivity
    e = zeros(Mesh.Jnodes,1); % emissivity
    rho_s_mass = zeros(Mesh.Jnodes,1);
    mprime = zeros(ssp,Mesh.Jnodes);
    p = zeros(Mesh.Jnodes,1); % non-staggered pressure points
    flux = zeros(Mesh.Jnodes,1); % staggered convective flux points
    yj = zeros(gsp,Mesh.Jnodes); % gas-phase species mass fraction
    j = zeros(gsp,Mesh.Jnodes); % diffusive flux
    D3 = zeros(Mesh.Jnodes,1); % diffusivity
    rhogphi = zeros(Mesh.Jnodes,1);
    drgpydt = zeros(gsp,Mesh.Jnodes);
    drhogphidt = zeros(Mesh.Jnodes,1);
    
    for i = 1:Mesh.Jnodes
        m(:,i) = yy1(ssp*(i-1)+1:ssp*(i-1)+ssp);
        m(:,i) = m(:,i)./MW;
    end
    
    yi = zeros(length(s_index),Mesh.Jnodes);
    T = yy1(ssp*Mesh.Jnodes+1:(ssp+1)*Mesh.Jnodes);
    rho_s_mass(:) = yy1((ssp+1)*Mesh.Jnodes+1:(ssp+2)*Mesh.Jnodes);
    rhogphi(:) = yy(1:Mesh.Jnodes);
    rgpy = reshape(yy(Mesh.Jnodes+1:end),gsp,Mesh.Jnodes);

    for i = 1:gsp
        yj(i,:) = rgpy(i,:)./trasspose(rhogphi(:));
    end
    
    %R = 8.314; 
    air = zeros(gsp,1); 
    air(end) = 1;
    
    for i = 1:Mesh.Jnodes
        yi(:,i) = m(s_index,i).*MW(s_index)./sum(m(s_index,i).*MW(s_index));
        phi(i) = phii(yi(:,i),rho_s_mass(i));
        k(:,i) = afac .*((T(i)).^nfac).* exp(-ea ./(R*T(i)));
        mprime(:,i) = reactions*(k(:,i).*m(istart,i)).*MW; % dm/dt
        wdot_mass(:,i) = mprime(:,i)./ Mesh.dv;
        kb(i) = kba(T(i),yi(:,i), phi(i),rho_s_mass(i)); % thermal conductivity [W/m/K]
        e(i) = epsilon(yi(:,i),rho_s_mass(i),phi(i));
        M = 1/sum(yj(:,i)./MW(g_index)); 
        p(i) = rhogphi(i)/phi(i)*R*abs(T(i))/M-p0;
        % diffusivity
        D3(i) = .018829*sqrt(T(i)^3*(1/32+1/28))/((p(i)+p0)*5.061^2*.93);
    end 
         
     for i = 1:Mesh.Jnodes-1
         flux(i) = -Kd*1/D3(i)*((p(i+1)-p(i))/Mesh.dz- rhogphi(i)/phi(i)*10*0);

        if flux(i)<0
            flux(i) = 0;
        end

        for k = 1:gsp
            D = D3(i+1)*D3(i)/(D3(i+1)+(D3(i)-D3(i+1))/2);
            rhophi = (rhogphi(i)+rhogphi(i+1))/2;
            j(k,i) = -D/MW(k)*rhophi*(yj(k,i+1)-yj(k,i))/Mesh.dz;
        end
     end
    
    flux(Mesh.Jnodes) = -Kd/D3(end)*((0-p(Mesh.Jnodes))/(Mesh.dz/2)...
        - rhogphi(Mesh.Jnodes)/phi(Mesh.Jnodes)*10*0);

    if flux(Mesh.Jnodes) < 0
        flux(Mesh.Jnodes) = 0;
    end

    for ii = 1:gsp
        j(ii,Mesh.Jnodes) = 0-.01*(air(ii)-yj(ii,end));
        if j(ii,Mesh.Jnodes) < 0
            j(ii,Mesh.Jnodes) = 0;
        end
    end
    
    for i = 2:Mesh.Jnodes-1
        drhogphidt(i) = sum(wdot_mass(g_index,i)) - (flux(i)-flux(i-1))/Mesh.dz;
        yfi = yj(:,i).*flux(i);
        yfii = yj(:,i-1).*flux(i-1);
        drgpydt(:,i) = wdot_mass(g_index,i) - (yfi-yfii)./Mesh.dz -(j(:,i)-j(:,i-1))/Mesh.dz;    
    end
    
    drhogphidt(1) = sum(wdot_mass(g_index,1)) - flux(1)/Mesh.dz;
    drgpydt(:,1) = wdot_mass(g_index,1) - (yj(:,1)*flux(1))./Mesh.dz - j(:,1)./Mesh.dz;
    drhogphidt(end) = sum(wdot_mass(g_index,end)) - (flux(end)-flux(end-1))/Mesh.dz;
    yfi = yj(:,end)*flux(end);
    yfii = yj(:,end-1)*flux(end-1);
    drgpydt(:,end) = wdot_mass(g_index,end) - (yfi-yfii)./Mesh.dz -(j(:,end)-j(:,end-1))/Mesh.dz;
    
    % mass flux at top boundary (where gases are released)
    tempflux = flux(Mesh.Jnodes)+sum(j(:,Mesh.Jnodes),'all');
    dydt = [drhogphidt(:); drgpydt(:)];  
end

% ODE function for energy trassport
function [dydt] = yprime1(t, yy, Mesh) %, reactions, afac, nfac, ea, istart, s_index, g_index, qs, MW, deltah, ssp)

    global reactions afac nfac ea istart s_index g_index qs MW deltah ssp 

    wdot_mass = zeros(ssp,Mesh.Jnodes);
    k = zeros(28,Mesh.Jnodes);
    m = zeros(ssp,Mesh.Jnodes);
    phi = zeros(Mesh.Jnodes,1);
    kb = zeros(Mesh.Jnodes,1);
    e = zeros(Mesh.Jnodes,1);
    rho_s_mass = zeros(Mesh.Jnodes,1);
    drhosdt = zeros(Mesh.Jnodes,1);
    mprime = zeros(ssp,Mesh.Jnodes);
    Tprime = zeros(Mesh.Jnodes,1);   
    
    for i = 1:Mesh.Jnodes
        m(:,i) = yy(ssp*(i-1)+1:ssp*(i-1)+ssp);
        m(:,i) = m(:,i)./MW;
    end
    
    yi = zeros(length(s_index),Mesh.Jnodes);
    T = yy(ssp*Mesh.Jnodes+1:(ssp+1)*Mesh.Jnodes);
    rho_s_mass(:) = yy((ssp+1)*Mesh.Jnodes+1:(ssp+2)*Mesh.Jnodes);
    % R = 8.314; % univ gas ct [J/molK]
    sigma = 5.670374419e-8; % sb const for rad [W/(m2K4)]
    h = 10; % heat transfer coefficient for heat xfr from heated surface to surroundings
    tr = 0;
    
    for i = 1:Mesh.Jnodes
        yi(:,i) = m(s_index,i).*MW(s_index)./sum(m(s_index,i).*MW(s_index));
        phi(i) = phii(yi(:,i),rho_s_mass(i));
        k(:,i) = afac .*((T(i)).^nfac).* exp(-ea ./(R*T(i)));
        mprime(:,i) = reactions*(k(:,i).*m(istart,i)).*MW; % dm/dt
        wdot_mass(:,i) = mprime(:,i)./ Mesh.dv;
        kb(i) = kba(T(i),yi(:,i), phi(i),rho_s_mass(i)); % thermal conductivity W/m/K
        e(i) = epsilon(yi(:,i),rho_s_mass(i),phi(i));
    end 
    
    for j = 2:Mesh.Jnodes-1
       ddd = cp(T(j));
       deltah = q_srxns(T(end));
       Tprime(j) = (1/(Mesh.dz^2)*((kb(j)+kb(j+1))/2*(T(j+1)-T(j))+(kb(j)+kb(j-1))/2*(T(j-1)-T(j)))...
           +e(j)*tr/Mesh.Jnodes*qs/Mesh.dz+sum(abs(wdot_mass(istart,j)).*q_srxns(T(j))))...
           /(rho_s_mass(j).*sum(ddd(s_index).*yi(:,j))); % dT/dt
    end

    de = cp(T(end));
    c = sum(de(s_index).*yi(:,end));
    Tprime(Mesh.Jnodes) = (Mesh.a*(e(end)*qs*(1-tr)-h*(T(end)-300)-e(end)*sigma*(T(end)^4-300^4))...
        -Mesh.a*(kb(end)+kb(end-1))/2*(T(end)-T(end-1))/Mesh.dz...
     +Mesh.dv*sum(abs(wdot_mass(istart,end)).*q_srxns(T(end))))/(Mesh.dv*rho_s_mass(end)*c);
    d1 = cp(T(1));
    Tprime(1) = (Mesh.a*kb(1)/(Mesh.dz)*(T(2)-T(1))+Mesh.a*e(1)*tr/Mesh.Jnodes*qs...
        +Mesh.dv*sum(abs(wdot_mass(istart,1)).*q_srxns(T(1))))...
        /(Mesh.dv*rho_s_mass(1)*sum(d1(s_index).*yi(:,1)));
     
    for i = 1:Mesh.Jnodes
        drhosdt(i) = - sum(wdot_mass(g_index,i));
    end

    dydt = [mprime(:); Tprime(:); drhosdt(:)];  
end

% conductivity function
function kb = kba(T, yi, phi, rho_s_mass, s_index, solid_densities)
    % [W/m/K]
    %global s_index MW

    k = zeros(length(s_index),1)+.17*(T/300)^.594;
    k(19) = .6;
    k(3) = .065*(T/300)^.435+5.670374419e-8*3.3e-3*(T)^3;
   
    % s_density.*MW(s_index)*1000;
    % = [9.37000000000000;9.37000000000000;25;11.5000000000000;11.5000000000000;...
    %     11.5000000000000;5.88000000000000;3.48000000000000;3.59000000000000;...
    %     5.88000000000000;4;7.29000000000000;5.76000000000000;7.22000000000000;5;1.67000000000000;...
    %     55;0.00369448575008421;0.00580475748165167;0.00541504238833635;0.0806563778419628;...
    %     0.0101350128633761;0.00507436386823035;0.00579578562602859; 1777777777].*MW(s_index)*1000;
    yi(18:end) = 0;
    kb = rho_s_mass*sum(yi.*k./(solid_densities))/(1-phi);
   
end

% emissivity function
function e = epsilon(yi, rho_s_mass, phi, s_index, solid_densities)

    %global s_index MW
    e = zeros(length(s_index),1)+0.757;
    e(3) = 0.957; % char
    e(19) = .95; % acqua
    e(49) = e(3); % ash

    
    % s_density.*MW(s_index)*1000;
    % = [9.37000000000000;9.37000000000000;25;11.5000000000000;11.5000000000000;...
    %     11.5000000000000;5.88000000000000;3.48000000000000;3.59000000000000;...
    %     5.88000000000000;4;7.29000000000000;5.76000000000000;7.22000000000000;5;1.67000000000000;...
    %     55;0.00369448575008421;0.00580475748165167;0.00541504238833635;0.0806563778419628;...
    %     0.0101350128633761;0.00507436386823035;0.00579578562602859; 177777777].*MW(s_index)*1000; 
    yi(18:end) = 0;
    e = rho_s_mass*sum(yi.*e./(solid_densities))/(1-phi);
end 

% heat capacity function
function cp = cp(ssp, T)
    %global ssp
    % cp in [J/kg/K]
    cp = zeros(ssp,1)+(1.5+.001*T)*1000;
    cp(15) = (.7+.0035*T)*1000;
    cp(39) = 4188; % water
    % cp(49) = ? for N2 %%%%%%%%%%%%%%%%%
    cp(50) = 880; % ASH
end

% heat of reactions function
function q_srxns = q_srxns(MW, istart, deltah, T)
    % J/kg of reactant
    %global MW istart
    
    %deltah = [-1300; 27100; 23200; -62700; -5000; -500; -42400; 17900; 12000; -10300; 30700; 26000; -31100;...
     %  -26100; 46200; -21100; -83600; 1300; 1300; 10100; -29100; -13400; 48600; 0; 0; 0; 0; 0; 10000; 10000; 10000]*4.184; % added reactions for char
    q_srxns = deltah*4.184./MW(istart); % convert to J/kmol
    q_srxns(28) = -2.41e6;
end

% porosity function
function phi = phii(yi, rho_s_mass, solid_densities)

    %global MW s_index
    
    %s_density = %[9.37000000000000;9.37000000000000;25;11.5000000000000;11.5000000000000;...
        %11.5000000000000;5.88000000000000;3.48000000000000;3.59000000000000;...
        %5.88000000000000;4;7.29000000000000;5.76000000000000;7.22000000000000;5;1.67000000000000;...
        %55;0.00369448575008421;0.00580475748165167;0.00541504238833635;0.0806563778419628;...
        %0.0101350128633761;0.00507436386823035;0.00579578562602859; 1777777]*1000; 
    % s_density*1000;    
    yi(18:end) = 0;
    phi = 1-sum(yi./(solid_densities))*rho_s_mass;
    % phi = .7432;
end