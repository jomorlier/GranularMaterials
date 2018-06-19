%% Impl?mentation de la methode de calcul transitoire "Verlet" , Prof Joseph Morlier, ISAE-SUPAERO
%
% Set system parameters
% E-k-m-k-m-k-m-xw
% We consider a linear chain of 3 particles of mass m = 1  connected to each other by a spring of stiffness k = 1. 
% The distance between the particles at the beginning is xe = 2. 
% The first and the last particle are connected to a wall: 
% the one on the left is fixed, the other oneis moved very slowly by a
% distance xw = xE.
% from its equilibrium position and then fixed in that position. 


%%
% soit $F(z)$ une fonction, et $x_0$ un r?el,
% 
function test_Nbody
%%J.Morlier
%Exemple Nbody fully parametrized for Optimization
% Limitation 1 dof par mass
clear all; close all;
grain=grain_init();

Domain=Domain_init(grain);

interaction=interaction_init(Domain);
 
core_init=core_init(grain,interaction,Domain);

[x,v,a,f]=core(grain,interaction,Domain,core_init);

[Ek,Ep,Etot]=Energies(v,x,grain,Domain,interaction,core_init);

Results=Results(x,v,a,f,Etot,Ek,Ep,core_init);

vizu(grain,Results,Domain);



function grain=grain_init()
%Set systegrain.m paragrain.meters
grain.id=1;
grain.m =1; %grain.mass
grain.E=1;
grain.R=0.1;

function Domain=Domain_init(grain)
Domain.Nm=10;%Number of Mass
Domain.connectivity= [1:Domain.Nm;2:Domain.Nm+1]';
Domain.BC.x0=0; %position of the walls
Domain.xe=2;
Domain.BC.xend=(Domain.Nm+2)*Domain.xe;
Domain.m= 1*ones(1,Domain.Nm) %1:Domain.Nm%1*ones(1,Domain.Nm);
Domain.mBC=[1 Domain.m 1];

function interaction=interaction_init(Domain)
interaction.Nk=Domain.Nm+2;
interaction.k=  1*ones(1,interaction.Nk)%1:interaction.Nk%1*ones(1,interaction.Nk); %stiffness spring
interaction.xe = 2;
interaction.mesh=1;
interaction.contact=0; %Hertz contact? displacement
interaction.c=0; % Damping proportionnal speed NL?

function core_init=core_init(grain,interaction,Domain)
%Set initial conditions
%x0 = Domain.BC.x0; xend = Domain.BC.xend; %position of the walls
core_init.x_zero=interaction.xe*(1:Domain.Nm);
%interaction.xe*(1:Domain.Nm);

core_init.v_zero=0*(1:Domain.Nm);
% impact=-10*grain.m;
% %% 
% core_init.v_zero(end)=impact;
%Set simulation parameters
core_init.tf =15; %time of simulation
core_init.delta_t =0.01; %time step
core_init.Nt = ceil(core_init.tf/core_init.delta_t); %number of time increments
core_init.t=(0:core_init.Nt)*core_init.delta_t;%time

function [x,v,a,f]=core(grain,interaction,Domain,core_init)
%Initialize the loop variables
for n = 1 : core_init.Nt+1
    if n == 1
        x(n,:) = [Domain.BC.x0, core_init.x_zero, Domain.BC.xend];
        v(n,:) = [0, core_init.v_zero, 0];
        for jj = 1 : interaction.Nk
            if jj == 1
                f(n,jj) = 0;
                a(n,jj) = f(n,:);
            elseif jj >=2 & jj<=interaction.Nk-1
                f(n,jj) = (interaction.k(jj) * (x(n,jj+1) - 2 * x(n,jj) + x(n,jj-1)));
               
            else
                f(n,jj) = 0;
                a(n,jj) = f(n,jj);
            end
        end
        
         a(n,:) = f(n,:)./ Domain.mBC;
        alfa(n,:) = x(n,:) + (core_init.delta_t * v(n,:)) + (0.5 * a(n,:) * ...
        (core_init.delta_t^2)); %uniformly accelerated motion
    else
        x(n,:) = alfa(n-1,:);
        for ii = 1 : interaction.Nk
            if ii == 1
                f(n,ii) = 0 ;
            elseif ii >=2 & ii<=interaction.Nk-1
            f(n,ii) =(interaction.k(ii) * (x(n,ii+1) - 2 * x(n,ii) + x(n,ii-1)));%Verlet 
          
            else
            f(n,ii) = 0;
            end
         end

a(n,:) = f(n,:)./ Domain.mBC;
alfa(n,:) = (2 * x(n,:)) - x(n-1,:) + ((core_init.delta_t^2) * a(n,:));
v(n,:) = (alfa(n,:) - x(n-1,:))/(2 * core_init.delta_t);
    end
end
%       %Analytical solution
% x_rel = x(:,2) - x(:,1) - interaction.xe; 
% x_an = -2:0.1:2;
% F1 = grain.k * x_an;

function [Ek,Ep,Etot]=Energies(v,x,grain,Domain,interaction,core_init)
%Conservation of energy

for n = 1 : core_init.Nt+1
    for mm=1:interaction.Nk-1
Ekk(mm,n,1) = 0.5 * Domain.mBC(mm) * (v(n,mm))^2; %grain.kinetic energy of grain.mass1 
Epp(mm,n,1)= (0.5 * interaction.k(mm) * (x(n,mm+1) - x(n,mm) - interaction.xe)^2);
    end

Ek(n,1)=sum(Ekk(:,n,1));
Ep(n,1)=sum(Epp(:,n,1));

 % system grain.kinetic energy 
Etot(n,1) = Ek(n,1) + Ep(n,1); %total energy
end



function Results=Results(x,v,a,f,Etot,Ek,Ep,core_init)
%save Results
Results.x=x;
Results.v=v;
Results.a=a;
Results.f=f;
Results.t=core_init.t;
Results.Etot=Etot;
Results.Ek=Ek;
Results.Ep=Ep;
save Results.grain.mat Results

function vizu(grain,Results,Domain)
%visu
x=Results.x(:,2:end);
v=Results.v(:,2:end);
a=Results.a(:,2:end);
f=Results.f(:,2:end);
t=Results.t;
Etot=Results.Etot;
Ek=Results.Ek;
Ep=Results.Ep;

figure;
plot(x,t);xlabel('position x');ylabel('time');

figure;
plot(t,v(:,2:end));ylabel('speed v');xlabel('time');

figure; 
plot(t,f(:,2:end));ylabel(' force f');xlabel('time');

figure;
plot(t,a(:,2:end));ylabel('acceleration a');xlabel('time');

figure;
plot(t,[Ek,Ep,Etot])
ylabel('Energies');xlabel('time');

%figure;
% subplot(1,2,1);hold on;
% h=line([0,10],[0,0],'color','r','marker', 'o', 'linewidth', 2);
% i = 1;
% while i<=length(x)
% plot(x,t);xlabel('x');ylabel('t');
% set(h, 'ydata', t(i)) % Moves the line to the time indicated by t
%  drawnow % necessary to get figure updated
% end
% subplot(1,2,2);hold on;
% Draw initial figure
% xlim([1,9]);
% ylim([-1.5,1.5]);

figure; 
set(gcf,'Renderer','OpenGL'); 
% Animation Loop
i = 1;

while i<=length(x)
for ip=1:Domain.Nm

h(ip) = plot(x(1,ip),0,'o','MarkerSize',20,'MarkerFaceColor','b');hold on;



    
    
    set(h(ip),'XData',x(i,ip)); 
    
  

    drawnow;
  
    i = i+1;
end

end



