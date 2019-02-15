function [c_aq, n_napl] = degrade(c_aq_old, n_napl_old, fcfb, fcfa, fbfa, poros, Si, rho, MW, dx, A)
%DEGRADE Oil degradation reaction
%  fa A + fb B -> fc C
% reaction of napl compounds with oxygen (electron donor) facilitated by
% microorganisms reaction

% initialize concentrations for equilibration
c_aq = c_aq_old;
n_napl = n_napl_old;

% initialize total concentrations for reaction
c_a_tot = c_aq_old(:,1:24);
c_b_tot = c_aq_old(:,25);

% initial aqueous concentrations for reaction 
c_a = c_aq_old(:,1:24);
c_b = c_aq_old(:,25);
c_c = zeros.*(c_aq_old(:,1:24));

% loop over all cells in model domain
for i =1:size(c_aq_old,1) 

    % skip over cell if oxygen or napl is depleted
    if c_b(i,:) <=0 || sum(c_a(i,:)) <= 0 
        % skip the iteration, on to next cell
        continue
    end
  
    it = 0;
    % enter reaction within cell if oxygen or napl is present
    while sum(c_a(i,:)) > 0 && c_b(i,:) > 0 &&  it < 100
        
    it = it + 1;
        
   % React hydrocarbons as ratio of compounds present in aqueous phase
    Xi = c_a(i,:)./sum(c_a(i,:));   
    
    % reaction that occurs if napl compounds are greater than oxygen    
    if sum(fcfa.*c_a_tot(i,:)) < sum(Xi.*(fcfb).*c_b_tot(i,:))
       c_a(i,:) = zeros.*c_a_tot(i,:);
       c_b(i,:) = c_b_tot(i,:) - sum(Xi.*fbfa.*c_a_tot(i,:));
       c_c(i,:) = Xi.*fcfa.*c_a_tot(i,:);
              
    % reaction that occurs if napl compounds are greater than oxygen
    else
       c_a(i,:) = c_a_tot(i,:) - Xi.*(1./fbfa).*c_b_tot(i,:);
       c_b(i,:) = zeros.*c_b_tot(i,:);
       c_c(i,:) = Xi.*fcfb.*c_b_tot(i,:);
    end
    
   c_a_tot(i,:) = c_a(i,1:24) + Xi.*(1./fcfa).*c_c(i,:);
   c_b_tot(i,:) = c_b(i,:) + sum(Xi.*(1./fcfb).*c_c(i,:));  
   
   end
   
   % equilibration following reaction 
    [a,n_napl(i,1:24)] = equilibrate(c_a(i,1:24), n_napl(i,1:24), poros, Si, rho, MW, dx, A);
    c_aq(i,1:24) = a;
    c_aq(i,25) = c_b(i,:);
    
end

end








