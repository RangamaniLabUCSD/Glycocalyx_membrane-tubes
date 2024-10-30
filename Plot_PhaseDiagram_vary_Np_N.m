clc;
clear;
close all;

%% Define parameters
%% membrane properties
kappa = 10;  %% bending rigidity, Unit: k_BT
sigma = 0.012;  %% membrane tension, Unit: k_BT/nm^2
lambda = 0.0;   %% Line tension, Unit: k_BT/nm
c_0 = 0.00;   %% spontaneous curveture, Unit: nm^(-1)
f = 0;
%% Polymer properties
a = 10;      %% monomer size,  Unit: nm
Chi = 0;   %% Flory parameter
v = (1-2*Chi)*a^3;     %% Excluded volume parameter, Unit: nm^3

R_0 = 100;        %%   membrane patch radius,  Unit: nm
A_s = pi*R_0^2;   %%   membrane patch area,  Unit: nm^2

% xi_t = [12 14.23 16 18 20];
xi_t = linspace(10,30,201);
% xi_t = 12;
t = 0;
u = 0;
for xi = xi_t
    xi
    u = u+1;
    S_R = xi^2;     %% local area per chain, Unit: nm^2
    rho_g = 1/xi^2;  %% grafting density, Unit: nm^(-2)
    N_p = A_s/S_R;   %% Number of polymer chain

    N_j = linspace(10,50,41);
    % N_j = [5 6 7 8 9 10];
    % N_j = 20;
    j = 0;
    for N = N_j
        N
        j = j+1;

        if xi<2*(1/sqrt(6))*a*N^(3/5)
            brush = 1;
        else
            brush = 0;
        end

        eta_span = 0.0:0.001:100; %% Shape parameter.
        i = 0;
        F_tot = 0;
        for eta = eta_span
            i = i+1;
            c_flat = (3/v/a^2/S_R^2)^(1/3);
            H_flat = N*S_R^(-1/3)*(v*a^2/3)^(1/3);
            F_flat = N*(9/2)*(v/3)^(2/3)*(S_R*a)^(-2/3)/pi/kappa;   %% energy of per chain on a flat surface
            R = R_0/(sqrt(2*(1+eta)));
            R_tube(i) = R_0/(sqrt(2*(1+eta)));
            L = R*eta;
            L_tube(i) = R*eta;
            H_tube = R*(1+4/3*H_flat/R)^(3/4)-R;
            thick_tube(i) = R*(1+4/3*H_flat/R)^(3/4)-R;
            N_ptube = 2*pi*R*L/xi^2;
            N_pcap = 2*pi*R^2/xi^2;
            F_poly1 = (9/2/kappa)*(2*R^3/xi^2)*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+5*N/3/R*(v*a^2/3/xi^2)^(1/3))^(1/5)-1)...
                +(9/2/kappa)*(eta*R^3/xi^2)*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+4*N/3/R*(v*a^2/3/xi^2)^(1/3))^(1/2)-1);
            F_poly2 = N_pcap*F_flat*R/H_flat*3*((1+5/3*H_flat/R)^(1/5)-1)+N_ptube*F_flat*R/H_flat*3/2*((1+4/3*H_flat/R)^(1/2)-1);
            F_line = 2*lambda/kappa*R;
            F(i) = (9/2/kappa)*(2*R^3/xi^2)*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+5*N/3/R*(v*a^2/3/xi^2)^(1/3))^(1/5)-1)...
                +(9/2/kappa)*(eta*R^3/xi^2)*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+4*N/3/R*(v*a^2/3/xi^2)^(1/3))^(1/2)-1)...
                + eta*R^2*(1/R-c_0)^2 + R^2*(2/R-c_0)^2 + sigma/kappa*R^2*(2*eta+1) + 2*lambda/kappa*R - f*R*(1+eta)/(pi*kappa);
            F_tot = F_tot+exp(-F(i));
        end
        %% Energy minimum
        [y,index_min] = min(F);
        eta_min(j) = eta_span(index_min);
        R_tube_min(j) = R_tube(index_min);
        L_tube_min(j) = L_tube(index_min);
        thick_tube_min(j) = thick_tube(index_min);
        E{j} = F;
        leg{j}=strcat('\xi=',num2str(xi),'nm');

        if eta_span(index_min)==0  
            phase = 0; %% No tubular
        else
            phase = 1; 
        end
        t = t+1;
        Results(t,:) = [xi,rho_g,N_p,N,eta_span(index_min),R_tube(index_min),L_tube(index_min),thick_tube(index_min),phase,kappa,sigma,lambda,c_0,f,a,R_0,brush];
    end
    %% find the boundary of each phase
    if eta_min(end)==0
        BC1 = find(eta_min==0);
        index_BC1 = BC1(end);
    else
        BC1 = find(eta_min>0);
        index_BC1 = BC1(1);
    end   
    BC(u,:) = [xi,rho_g,N_p,N_j(index_BC1)];
end
eta_span = eta_span';
subplot(2,2,1);
for j=1:length(N_j)
    plot(eta_span,E{j},'linewidth',2)
    hold on
end
hold off
legend(leg)
xlabel('\eta');
ylabel('\itF_{tot}/\pi\kappa')
str = {'\xi=15 nm','\lambda=0.5 k_BT/nm','c_0=0.04 nm^{-1}','a=10 nm','N=20','R_0=50 nm'};
subplot(2,2,2);
plot(N_j,eta_min,'r','linewidth',2)
xlabel('\itN');
ylabel('\it\eta_{min}')
subplot(2,2,3);
scatter(Results(:,3),Results(:,4),[],Results(:,5),'filled')
u=colorbar;
set(u,'FontName','Times New Roman','FontSize',15,'linewidth',2.5,'FontWeight','bold');
set(get(u,'title'),'string','\eta');
% u.Label.String = '\eta';
set(gca,'xscale','log');
xlabel('\itN_p');
ylabel('\itN')
box on
% plot(N_j,R_tube_min,'b','linewidth',2)
% xlabel('\itN');
% ylabel('\itR_{min}')
subplot(2,2,4);
plot(N_j,L_tube_min,'g','linewidth',2)
xlabel('\itN');
ylabel('\itL_{min}')