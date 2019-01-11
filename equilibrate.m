function [c_aq,n_napl] = equilibrate(c_aq_old, n_napl_old, poros, Si, rho,...
    MW, dx,A)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% initialize with old values
c_aq = c_aq_old;
n_napl= n_napl_old;

for i=1:size(c_aq_old,1)
    % evaluate total moles of compounds in cell
    n_tot = n_napl(i,:) + c_aq_old(i,:) * (poros*dx*A-sum(n_napl_old(i,:).*MW./rho));
    
    % make sure that while loop is entered
    n_napl_temp = Inf(size(n_napl(i,:)));
    
    % Picard iteration until nothing changes anymore
    it = 1;
    while (sum(sqrt((n_napl(i,:) - n_napl_temp).^2))  >    1e-6) && it < 100
        it = it+1;
        % save step for convergence criterion
        n_napl_temp = n_napl(i,:);
        % Picard equation
        n_napl(i,:) = n_tot./(1+Si./sum(n_napl(i,:)) * (poros*dx*A-sum(n_napl(i,:).*MW./rho)) ) ;
    end
    if it > 100
        warning('Iteration exceeded 100 calculations.');
    end
    % get aqueous concentrations
    c_aq(i,:) = Si.*n_napl(i,:)./sum(n_napl(i,:));
end
end