% Plotting of breakthrough curves
% Plotting at different time points using data 
% Individual subplots for hydrocarbon type 


loop_y = [0.01 0.1 0.2 0.3 0.4 0.5 1 2 5 10 15 20 25 30 35 40];  % Enter time in years you would like to plot
loop_s = (loop_y.*365.*24.*3600);           % Conversion to seconds

for i = 1:length(loop)
% Alkane subplot 
set(groot,'defaultAxesColorOrder',co([1:8,25],:))
figure(2)
subplot(2,2,1)
plot([dt:dt:loop_s(:,i)]./(3600.*24.*365),BTC([1:loop_s(:,i)./dt],[1:8,25])./[c_aq_i(1,1:8),c_in(:,25)]);
xlabel('t [years]')
ylabel('c/c_0 [-]')
% legend(compound(:,[1:8,25]), 'Position', [0.3507 0.6528 0.1103 0.2030])
legend(compound(:,[1:8,25]), 'Location', 'eastoutside')
title('Breakthrough Curve - Alkanes')

% Alkene subplot
set(groot,'defaultAxesColorOrder',co([9:15,25],:))
subplot(2,2,2)
plot([dt:dt:loop_s(:,i)]/(3600.*24.*365),BTC([1:loop_s(:,i)./dt],[9:15,25])./[c_aq_i(1,9:15),c_in(:,25)])
xlabel('t [years]')
ylabel('c/c_0 [-]')
% legend(compound(:,[9:15,25]), 'Position', [0.8387 0.6638 0.1103 0.1810])
legend(compound(:,[9:15,25]), 'Location', 'eastoutside')
title('Breakthrough Curve - Alkenes')

% BTEX subplot 
set(groot,'defaultAxesColorOrder',co([16:19,25],:))
subplot(2,2,3)
plot([dt:dt:loop_s(:,i)]/(3600.*24.*365),BTC([1:loop_s(:,i)./dt],[16:19,25])./[c_aq_i(1,16:19),c_in(:,25)])
xlabel('t [years]')
ylabel('c/c_0 [-]')
% legend(compound(:,[16:19,25]), 'Position', [0.3804 0.2230 0.1103 0.1149])
legend(compound(:,[16:19,25]), 'Location', 'eastoutside')
title('Breakthrough Curve - BTEX')

% Additive subplot 
set(groot,'defaultAxesColorOrder',co([20:25],:))
subplot(2,2,4)
plot([dt:dt:loop_s(:,i)]/(3600.*24.*365),BTC([1:loop_s(:,i)./dt],20:25)./[c_aq_i(1,20:24),c_in(:,25)])
xlabel('t [years]')
ylabel('c/c_0 [-]')
% legend(compound(:,20:25), 'Position', [0.9 0.10 0.1103 0.1369])
legend(compound(:,20:25), 'Location', 'eastoutside')
title('Breakthrough Curve - Additives')
save('data.mat');


filename=fullfile(DirectoryPath,['Breakthrough Curve_', num2str(loop_y(:,i)),' Years.jpeg']);
saveas(gcf, filename);
end
