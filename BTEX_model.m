% =========================================================================
% 1-D Simulation of Transport and Fate of BTEX 
% Hannah McNeil and Nicole Smith
% University of Tübingen
%
% December, 2018 
% =========================================================================

clear all

% Transport Coefficients
m_tot = 8000;       % total mass of NAPL [kg]
c_in = 0.008;       % Dissolved oxygen inflow concentration [kg/m3]

L = 10;             % length of aquifer [m] 
W = 20;             % width of aquifer [m]
H = 2;              % depth [m] 
poros = 0.35;       % porosity [-]
q = 1.1574e-05;     % Specific Discharge [m/s]
alpha = 0.01;       % Longitudinal Dispersivity [m]
Dp = 1e-9;          % pore diffusion coefficient [m2/s] -- assumed

tin = 0;            % begin time [s]
te = 3600*(1*240);  % end time [s]

% Derived Coefficients
A = H*W;            % area [m2]
v = q/poros;        % seepage velocity [m/s]
D = alpha*v+Dp;     % dispersion coefficient [m2/s]

% Individual gasoline compounds and oxygen 
compound = {'Propane' 'i-Butane / n-Butane' 'i-Pentane / n-Pentane' ...
    'i-Hexane / n-Hexane' 'i-Heptane / n-Heptane' 'i-Octane / n-Octane' ...
    'Nonane' 'Decane' 'Butene' 'Pentene' 'Hexene' 'Heptene' 'Octene' ...
    'Nonene' 'Decene' 'Benzene' 'Toluene' 'o-, p-Xylene' 'Ethylbenzene' ...
    'Trimethylbenzene' 'MTBE' 'Ethanol' 'Aniline' 'Phenol' 'Oxygen'};

% Individual gasoline compounds - chemical parameters

% density [kg/m3]
rho = [493 573 626 660 684 703 718 730 588 626 673 697 715 729 741 874 ...
    865 868 867 873 776 780 1022 1070];

% molar mass [kg/mol]
MW = [0.04409 0.05814 0.07215 0.08618 0.1002 0.11423 0.12825 0.14828 ...
    0.0561 0.07013 0.08417 0.09819 0.11222 0.12624 0.14019 0.07811 ...
    0.09214 0.10617 0.1062 0.12019 0.08815 0.0461 0.09313 0.09411];

% solubility constants at 20-25 C [kg/m3] 
Si = [0 0.061 0.040 0.013 0.0022 0.00066 0.0002 9.0e-06 0.221 0.148 ...
    0.05 0.0182 0.0041 0.00112 0.000115 1.79 0.526 0.162 0.15 0.0752 ...
    42 1000 36 82.80];

% Invidual gasoline compounds - initial conditions

% weight percent of each compound in initial mixture [fraction]
wt = [0 0.047 0.178 0.173 0.068 0.079 0.013 0.025 0.025 ... 
    0.0330 0.022 0.009 0.004 0.012 0.012 0.007 0.069 0.055 0.014 ...
    0.093 0.011 0.051 8.0e-06 3.7e-05];

% Invidual gasoline compounds - derived initial conditions 

% moles of each compound in initial mixture [mol]
n = wt*m_tot./MW;

% Xo = (n_i/n_total moles) molar fraction for mixture [-]
Xo = n./sum(n);

% c_aq = (Xo*Si) initial water concentration at equilibrium  [kg/m3]
c_aq_i = Xo.*Si;

% Requirement 1: Neumann-Number = 1/4
% Requirement 2: Courant-Number Cr = dt*v/dx = 1
% Solve for dx:
dx = 4*D/v;
% Solve for dt
dt = dx/v;

x=0.5*dx:dx:L;
nx = length(x);

% Number of Components [given compounds + oxygen] 
ncomp = 25;

% Matrix of Aqueous Concentrations - Compound and Oxygen
% Rows related to length coordinates
% Columns related to components
%c_aq = zeros(nx,ncomp);
c_aq = ones(nx,ncomp).*[c_aq_i,c_in(:,end)];

% Matrix of Initial Total NAPL Mass - Compound and Oxygen
% Rows related to length coordinates
% Columns related to components
m = ones(nx,1)*[(m_tot.*dx/L.*wt), 0];

% Matrix of Initial free phase NAPL Mass - Compound and Oxygen
% Rows related to length coordinates
% Columns related to components
%m_napl = zeros(nx,ncomp);

% Matrix of Inflow Concentrations - Compound and Oxygen 
% Columns related to components
c_in = [zeros(1,ncomp-1), c_in];

% Initial estimate of area for water flow 
%Aw = ones(nx,ncomp).*poros.*H.*W;

% initialize breakthrough curve
BTC=zeros(0,ncomp);

% Open figure and delete its content
figure(1);clf
% 
% % Open video 
% v = VideoWriter('transport_model.avi');
% open(v);

% =========================================================================
% Loop over all timepoints
% =========================================================================
for t=dt:dt:te

    BTC=[BTC;c_aq(end,:)];
%     
%     % =====================================================================
%     % EQUILIBRATION
%     % =====================================================================
%     [m_napl] = equilibrate3(m(:,1:ncomp-1),c_aq(:,1:ncomp-1));
% %     m_napl(m_napl<0) = 0;
%     
%     % update aqueous concentration and mass 
%     Aw = (A.*dx.*poros - ones(nx,ncomp).*sum((m_napl./rho),2))/dx;
%     c_aq(:,1:end-1) = ones(nx,ncomp-1).*Si.*m_napl(:,1:ncomp-1)./...
%         (ones(nx,1).*MW).*(ones(1,ncomp-1).*(sum(m_napl(:,1:ncomp-1)./...
%         ones(nx,1).*MW,2)));
    
    % =====================================================================
    % ADVECTION
    % =====================================================================
    % Advection a Courant-number 1 implies that the concentrations are
    % moved by exactly one box. The values in the last box are moved out.
    % The first box receives the inflow concentration.
    
    c_aq(2:end,:)= c_aq(2:end,:) + dt*(q/dx)*(c_aq(1:end-1,:) - c_aq(2:end,:));
    c_aq(1,:)= c_aq(1,:) + dt*(q/dx)*(c_in - c_aq(1,:));
    
    % =====================================================================
    % DISPERSION
    % =====================================================================

    Jd=(c_aq(1:end-1,:)-c_aq(2:end,:))/dx*D;
    % add a dispersive flux of zero at the inflow boundary and assume that
    % the dispersive flux at the outflow is identical to that at the last
    % internal interface
    Jd =[zeros(1,ncomp);Jd;Jd(end,:)];
    
    % concentration change due to divergence of dispersive flux 
    c_aq = c_aq + dt/dx*(Jd(1:end-1,:)-Jd(2:end,:));
    
    % =====================================================================
    % REACTION
    % =====================================================================
%    % =====================================================================
%     % EQUILIBRATION
%     % =====================================================================
%     [m_napl] = equilibrate2(m(:,1:ncomp-1),c_aq(:,1:ncomp-1));
%     m_napl(m_napl<0) = 0;
%     
%     % update aqueous concentration and mass 
%     Aw = A.*dx.*poros - ones(nx,ncomp).*sum((m_napl./rho),2);
%     c_aq(:,1:end-1) = ones(nx,ncomp-1).*Si.*m_napl(:,1:ncomp-1)./...
%         (ones(nx,1).*MW).*(ones(1,ncomp-1).*(sum(m_napl(:,1:ncomp-1)./...
%         ones(nx,1).*MW,2)));
%    
    % =====================================================================
    % GRAPHICAL OUTPUT
    % =====================================================================
   
    % Graphical output every 10 minutes
    figure(1)
    semilogy(x,c_aq);
    xlabel('x [m]');
    ylabel('c [g/L]');
    ylim([0 500]);
    legend(compound)
    title(sprintf('Concentration, t=%6.1fh',t/3600));
    drawnow
    
%     % Graphical video 
%         frame = getframe(gcf);
%         writeVideo(v,frame)
    
end

% % End video 
% close(v) 

figure(3)
%subplot(2,1,2)
plot([dt:dt:te]/3600,BTC)
xlabel('t [h]')
ylabel('c [g/L]')
ylim([0 100]);
title(sprintf('Breakthrough Curve'))% n = %3.1f')),n))
saveas(gcf, 'Breakthrough Curve.jpeg')





