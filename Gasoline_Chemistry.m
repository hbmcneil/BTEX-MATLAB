% Hannah McNeil and Nicole Smith 
% Transport and Fate of BTEX in an aquifer 
% Task One 
% November 2018 

% ======================================================================= %
% STEP ONE: Input all componenets and variables % 
% ======================================================================= %

% Aquifer Variables 
m_tot = 8000;       % total mass of NAPL [kg]

% Gasoline Component Variables 
% rho = density [kg/m3]
% MW = molar mass [kg/mol]
% Si = solubilities in mol/m3
% wt = weight percent in napl [-]
[compound, rho, MW, Si, wt] = BTEX_data();

% Carbon content for each compound by from the molecular formula 
C = [3; 4; 5; 6; 7; 8; 9; 10; 4; 5; 6; 7; 8; 9; 10; 6; 7; 8; 8; 9; 5; 2;
    6; 6];
C = C';

% Hydrogen content for each compound by from the molecular formula 
H = [8; 10; 12; 14; 16; 18; 20; 22; 8; 10; 12; 14; 16; 18; 20; 6; 8; 10;
    10; 12; 12; 6; 7; 6];
H = H';

% Oxygen content for each compound by from the molecular formula 
O = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 1; 0; 
    1];
O = O';

% Oxygen content for each compound by from the molecular formula 
N = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 
    0];
N = N';

% ======================================================================= %
% STEP TWO: Calculating new variables % 
% ======================================================================= %

n = wt*m_tot./MW;   % moles of each compound in initial mixture [mol] 
Xo = n./sum(n);     % molar fraction [%]

Cw_0 = Xo.*Si;      % Calculation of maximum solubility [mol/m3]

% Oxidizing Conditions
% Calculation of chemical oxygen demand 

% Calculation of COD proportion to individual components
COD_prop = C + (H./4) - (O./2) + (5.*N./4);

% Calculation of COD mass using initial mass (g)
COD_mass = (((wt./100).*m_tot.*1000)./MW).*COD_prop.*15.994.*2;

% Reducing Conditions
% Calculation of mineralogical Fe(III) demand

% Calculation of Fe proportion to individual components
Fe_prop = (C.*4) + H - (O.*2) + (5.*N);

% Calculation of COD mass using initial mass (g)
Fe_mass = (((wt./100).*m_tot.*1000)./MW).*Fe_prop.*55.845;


% print results and finish calculation for total FE required and time 

% ======================================================================= %
% Task 4: Calculating Plume Length
% ======================================================================= %

% do I need to use seperate variables or the same as btexmodel? 

% M = H or z of aquifer [m]
M = 2; 

%vertical transverse dispersivity [L] -- pg.12 EM2 5-20x smaller than alpha
%[chose 10x smaller]
alpha_t = 0.001;

%concentration of donor = NAPL concentration [mg/L] **check if mg/L is unit
%or mol/m3
c_d0 = c_aq(:,1:24);

%concentration of acceptor = oxygen concentration [mg/L]
c_a0 = 8;

%stoichiometric ratio [-]
gamma = C./COD_prop;

% Length of Plume if pi^2*alpha_l*alpha_t/M^2 << 1
if pi^2*alpha*alpha_t/M^2 < 1
    L_plume = ((4*M^2)/(pi*alpha_t))*(log((4/pi)*((gamma.*c_d0 + c_a0)/c_a0)));
end




