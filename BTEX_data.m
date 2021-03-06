function [compound, rho, MW, Si, wt, carbon, fcfb, fcfa, fbfa] = BTEX_data()
%BTEX_data Summary of gasoline data
% Chemical parameters for data calculations 

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

% solubilities in mol/m3
Si = Si./MW;

% weight percent of each compound in initial mixture [fraction]
wt = [0 0.047 0.178 0.173 0.068 0.079 0.013 0.025 0.025 ... 
    0.0330 0.022 0.009 0.004 0.012 0.012 0.007 0.069 0.055 0.014 ...
    0.093 0.011 0.051 8.0e-06 3.7e-05];

% Carbon content for each compound by from the molecular formula 
carbon = [3; 4; 5; 6; 7; 8; 9; 10; 4; 5; 6; 7; 8; 9; 10; 6; 7; 8; 8; 9; 5; 2;
    6; 6];
carbon = carbon';

% Hydrogen content for each compound by from the molecular formula 
hydrogen = [8; 10; 12; 14; 16; 18; 20; 22; 8; 10; 12; 14; 16; 18; 20; 6; 8; 10;
    10; 12; 12; 6; 7; 6];
hydrogen = hydrogen';

% Oxygen content for each compound by from the molecular formula 
oxygen = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 1; 0; 
    1];
oxygen = oxygen';

% Nitrogen content for each compound by from the molecular formula 
nitrogen = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 
    0];
nitrogen = nitrogen';

% Calculation of COD proportion to individual components
COD_prop = carbon + (hydrogen./4) - (oxygen./2) + (5.*nitrogen./4);

% Stoichiometric ratio between CO2 and O2
fcfb = carbon./COD_prop;

% Stoichiometric ratio between CO2 and napl compound 
fcfa = carbon;

fbfa = COD_prop;


end

