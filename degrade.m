function [c_aq, n_napl] = degrade(c_aq_old, n_napl_old, fcfb, fcfa, fbfa, poros, Si, rho, MW, dx, A)
%DEGRADE Oil degradation reaction
%  fa A + fb B -> fc C
% reaction of napl compounds with oxygen (electron donor) facilitated by
% microorganisms reaction

% initialize concentrations
c_aq = c_aq_old;
n_napl = n_napl_old;

c_a_tot = c_aq_old(:,1:24);
c_b_tot = c_aq_old(:,25);
c_a = c_a_tot;
c_b = c_b_tot;
c_c = zeros.*(c_aq_old(:,1:24));

% loop over all cells in model domain
for i =1:size(c_aq_old,1) 
    % skip over cell if oxygen or napl is depleted
     Xi = c_a_tot(i,:)./sum(c_a_tot(i,:));
    
    if c_b_tot(i,:) <=0 || sum(c_a_tot(i,:)) <= 0 
        % skip the iteration, on to next cell
        continue
    end
        
    % enter reaction within cell if oxygen or napl is present
   while sum(c_a(i,:)) > 0 && c_b(i,:) > 0 
    
    if sum(fcfa.*c_a_tot(i,:)) < sum(Xi.*fbfa.*c_b_tot(i,:))
       c_a(i,:) = zeros.*c_a_tot(i,:);
       c_b(i,:) = sum(Xi.*c_b_tot(i,:) - fbfa.*c_a_tot(i,:));
       c_c(i,:) = Xi.*fcfa.*c_a_tot(i,:);
       
    % enter reaction for second scenario within cell if oxygen or napl is present
    else
       c_a(i,:) = c_a_tot(i,:) - Xi.*(1./fbfa).*c_b_tot(i,:);
       c_b(i,:) = zeros.*c_b_tot(i,:);
       c_c(i,:) = Xi.*fcfb.*c_b_tot(i,:);
    end
    
    c_a_tot(i,:) = c_a(i,:) + (1./fcfa).*c_c(i,:);
    c_b_tot(i,:) = c_b(i,:) + sum((1./fcfa).*c_c(i,:));
    
   end
   
   [a,n_napl] = equilibrate(c_a(i,1:24), n_napl, poros, Si, rho, MW, dx, A);
    c_aq(i,1:24) = a;
    c_aq(i,25) = c_b(i,:);
end

end






