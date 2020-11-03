function [p_fit] = auto_fit(iso, parameters, run_op, spec_info, glacial_start)

%% Auto fit function
% Use the dD and d18o spectra structures with the nonlinear least squares
% fit routine to match a Gaussian (or summation of Gaussians) to each
% spectra window

% Assign depths/ages of each window
age = iso.P(1,:);  
% half-cm sample size
ss = .005;     

if run_op.q_1gauss == 1
p_fit = zeros(size(iso.P,2),5,3);    % initialize matrix to hold fit values
elseif run_op.q_2gauss == 1
p_fit = zeros(size(iso.P,2),7,3);    % initialize matrix to hold fit values
end
p_fit(:,1,1) = age;                 % plug in depth values

iso_mat3d = iso.P;                  % create 3d matrix with Power and confidence levels
iso_mat3d(:,:,2) = iso.cu;          % to allow fit for 95% confidence of power estimate
iso_mat3d(:,:,3) = iso.cl;          % as well

% calculate uncertainty of sigma if set to do so
if run_op.uncertainty == 1
    m_range = 1:3;
elseif run_op.uncertainty == 0
    m_range = 1;
end

% calculate fit for each spectra
for m = m_range
for n = 1:size(iso.P,2)
    
    % if spectra contains nans (ie is d17o before there is data), do not
    % attempt to fit this spectra
    if sum(isnan(iso.P(:,n))) > 0
        continue
    end
    
k = 2*pi*iso.f(3:end);     % import frequency data
f = k/(2*pi);               % define frequency vector

    if m == 1
        G_data = iso_mat3d(3:end,n,m);   % import power data, ignore depth and first 2 power values
    % set up ranges for uncertainty
    elseif m == 2      
        % for holocene spectra
        if age(n) <= glacial_start
        divide_start_up = find(f>=5,1,'first');
        divide_end_up = find(f>=7,1,'first');
        G_data = [iso_mat3d(3:divide_start_up-1,n,2);...
            linspace(iso_mat3d(divide_start_up,n,2),iso_mat3d(divide_end_up,n,3),divide_end_up-divide_start_up)';...
            iso_mat3d(divide_end_up:end,n,3)];
        % for glacial spectra
        elseif age(n) >glacial_start
        divide_start_up = find(f>=8,1,'first');
        divide_end_up = find(f>=11,1,'first');
        G_data = [iso_mat3d(3:divide_start_up-1,n,2);...
            linspace(iso_mat3d(divide_start_up,n,2),iso_mat3d(divide_end_up,n,3),divide_end_up-divide_start_up)';...
            iso_mat3d(divide_end_up:end,n,3)];
        end
    elseif m == 3
       % holocene
        if age(n) <= glacial_start
        divide_start_up = find(f>=5,1,'first');
        divide_end_up = find(f>=7,1,'first');
        G_data = [iso_mat3d(3:divide_start_up-1,n,3);...
            linspace(iso_mat3d(divide_start_up,n,3),iso_mat3d(divide_end_up,n,2),divide_end_up-divide_start_up)';...
            iso_mat3d(divide_end_up:end,n,2)];
        % glacial
        elseif age(n) >glacial_start
        divide_start_up = find(f>=8,1,'first');
        divide_end_up = find(f>=11,1,'first');
        G_data = [iso_mat3d(3:divide_start_up-1,n,3);...
            linspace(iso_mat3d(divide_start_up,n,3),iso_mat3d(divide_end_up,n,2),divide_end_up-divide_start_up)';...
            iso_mat3d(divide_end_up:end,n,2)];
        end
    end
   
    


% % even log spacing code from Tyler Jones
% % wavelength band to test - full range of frequencies
% freq_high = 1/100;     %meters
% freq_low = 1/2;   %meters
% 
% % number of evenly spaced logarithmic points to test withn wavelength band
% % this can cause an error if too large
% n_log = 16;
% 
% % equal log-frequency spacing, for power law fit
%  
% % determine b, for 10^b, for each wavelegth
%     b_high = log10(abs(1/freq_low));
%     b_low = log10(abs(1/freq_high));
%  
%     % Create a vector of 'n_log' logarithmically spaced points in the interval [10^b_high,10^b_low].
%     logspacing = logspace(b_high,b_low,n_log);
%     
%     % preallocate vectors
%     f_eq = ones(n_log-1,1);
%     P_eq = ones(n_log-1,1);
%     
%     for i=1:n_log-1
%         i_eq = find(f >= logspacing(i) & f < logspacing(i+1));    
%         f_eq(i) = mean(f(i_eq));
%         P_eq(i) = mean(G_data(i_eq));
%     end
% % end even log spacing code from Tyler Jones

% temporary assignment to debug the rest of code
f_eq = f;
k_eq = 2*pi*f_eq;   % logspace wavenumber
P_eq = G_data;
    
% remove annual peak if specified
if run_op.remove_peak ==1
    % remove annual peak - NEED TO MAKE THIS WORK FOR DEPTH WORLD
    % cut out one-year signal to find more accurate regression
    f_conc = f_eq;
    peak_start = nan(size(iso.P,2),1);
    peak_end = nan(size(iso.P,2),1);
    peak_start(1) = find (f_conc >= 3.8,1,'first'); % find beginning of annual peak
    peak_end(1) = find (f_conc >= 4.7,1,'first'); % find end of annual peak
    peak_start(2) = find (f_conc >= 3.6,1,'first'); 
    peak_end(2) = find (f_conc >= 4.7,1,'first');
    peak_start(3) = find (f_conc >= 3.6,1,'first'); 
    peak_end(3) = find (f_conc >= 5.05,1,'first');
    peak_start(4) = find (f_conc >= 3.6,1,'first'); 
    peak_end(4) = find (f_conc >= 5.05,1,'first');
    peak_start(5) = find (f_conc >= 3.6,1,'first'); 
    peak_end(5) = find (f_conc >= 5.1,1,'first');
    peak_start(6) = find (f_conc >= 3.6,1,'first'); 
    peak_end(6) = find (f_conc >= 5.45,1,'first');
    peak_start(7) = find (f_conc >= 3.6,1,'first'); 
    peak_end(7) = find (f_conc >= 5.4,1,'first');
    peak_start(8) = find (f_conc >= 3.6,1,'first'); 
    peak_end(8) = find (f_conc >= 5.7,1,'first');
    peak_start(9) = find (f_conc >= 4.7,1,'first'); 
    peak_end(9) = find (f_conc >= 6.4,1,'first');
    peak_start(10) = find (f_conc >= 4.7,1,'first'); 
    peak_end(10) = find (f_conc >= 6.9,1,'first');
    peak_start(11) = find (f_conc >= 4.9,1,'first'); 
    peak_end(11) = find (f_conc >= 6.7,1,'first');
    peak_start(12) = find (f_conc >= 4.9,1,'first'); 
    peak_end(12) = find (f_conc >= 5.0,1,'first');
    peak_start(13) = find (f_conc >= 5.8,1,'first'); 
    peak_end(13) = find (f_conc >= 7.3,1,'first'); 
    peak_start(14) = find (f_conc >= 5.9,1,'first'); 
    peak_end(14) = find (f_conc >= 7.6,1,'first');
    peak_start(15) = find (f_conc >= 7,1,'first'); 
    peak_end(15) = find (f_conc >= 8.5,1,'first');
    peak_start(16) = find (f_conc >= 7.4,1,'first'); 
    peak_end(16) = find (f_conc >= 9.1,1,'first');
    peak_start(17) = find (f_conc >= 8.2,1,'first'); 
    peak_end(17) = find (f_conc >= 9.5,1,'first');
    peak_start(18) = find (f_conc >= 4.9,1,'first'); 
    peak_end(18) = find (f_conc >= 5.0,1,'first');
    peak_start(19) = find (f_conc >= 4.9,1,'first'); 
    peak_end(19) = find (f_conc >= 5.0,1,'first');
    peak_start(20) = find (f_conc >= 4.9,1,'first'); 
    peak_end(20) = find (f_conc >= 5.0,1,'first');
    peak_start(21) = find (f_conc >= 4.9,1,'first'); 
    peak_end(21) = find (f_conc >= 5.0,1,'first');
    peak_start(22) = find (f_conc >= 8.9,1,'first'); 
    peak_end(22) = find (f_conc >= 10.7,1,'first');
    peak_start(23) = find (f_conc >= 7.2,1,'first'); 
    peak_end(23) = find (f_conc >= 9.8,1,'first');
    peak_start(24) = find (f_conc >= 4.9,1,'first'); 
    peak_end(24) = find (f_conc >= 5.0,1,'first');
    peak_start(25) = find (f_conc >= 4.9,1,'first'); 
    peak_end(25) = find (f_conc >= 5.0,1,'first');
    peak_start(26) = find (f_conc >= 9.2,1,'first'); 
    peak_end(26) = find (f_conc >= 11.2,1,'first');
    peak_start(27:end) = ones(1,size(iso.P,2)-26)*find(f_conc >= 4.9,1,'first'); 
    peak_end(27:end) = ones(1,size(iso.P,2)-26)*find(f_conc >= 5.0,1,'first');
%     peak_start(28) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(28) = find (f_conc >= 5.0,1,'first');
%     peak_start(29) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(29) = find (f_conc >= 5.0,1,'first');
%     peak_start(30) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(30) = find (f_conc >= 5.0,1,'first');
%     peak_start(31) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(31) = find (f_conc >= 5.0,1,'first');
%     peak_start(32) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(32) = find (f_conc >= 5.0,1,'first');
%     peak_start(33) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(33) = find (f_conc >= 5.0,1,'first');
%     peak_start(34) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(34) = find (f_conc >= 5.0,1,'first');
%     peak_start(35) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(35) = find (f_conc >= 5.0,1,'first');
%     peak_start(36) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(36) = find (f_conc >= 5.0,1,'first');
%     peak_start(37) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(37) = find (f_conc >= 5.0,1,'first');
%     peak_start(38) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(38) = find (f_conc >= 5.0,1,'first');
%     peak_start(39) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(39) = find (f_conc >= 5.0,1,'first');
%     peak_start(40) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(40) = find (f_conc >= 5.0,1,'first');
%     peak_start(41) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(41) = find (f_conc >= 5.0,1,'first');
%     peak_start(42) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(42) = find (f_conc >= 5.0,1,'first');
%     peak_start(43) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(43) = find (f_conc >= 5.0,1,'first');
%     peak_start(44) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(44) = find (f_conc >= 5.0,1,'first');
%     peak_start(45) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(45) = find (f_conc >= 5.0,1,'first');
%     peak_start(46) = find (f_conc >= 4.9,1,'first'); 
%     peak_end(46) = find (f_conc >= 5.0,1,'first');
%     f_conc(peak_start(n):peak_end(n)) = []; % remove annual peak
%     data_conc = G_data; % rename ln(P) vector
%     data_conc(peak_start(n):peak_end(n)) = []; % remove annual peak
    peak_length = length(f_conc(peak_start(n):peak_end(n))); % data points in peak
    P_eq(peak_start(n):peak_end(n)) = linspace(P_eq(peak_start(n)),P_eq(peak_end(n)),peak_length) ; % remove annual peak
end

% Single Gaussian fit    
if run_op.q_1gauss == 1
    
% define model for fit
G = @(p) p(1)*exp(-k_eq.^2*p(2)^2)...      % Gaussian diffusion function
    + p(3).^2*ss./abs(1-p(4)*exp(-1i*k_eq*ss)).^2; % Red noise function
    
% calculate residual
residual = @(p) log10(G(p)) - log10(P_eq);

%Initial parameters and limits for 1 gaussian fit
% for holocene spectra
if age(n) <= glacial_start
p_i = parameters(1,1:4,1);     % [P_1, sigma_1, sigma_noise, ar1]
p_lb = parameters(2,1:4,1);
p_ub = parameters(3,1:4,1);
% for glacial spectra
elseif age(n) > glacial_start
p_i = parameters(1,1:4,2);     % [P_1, sigma_1, sigma_noise, ar1]
p_lb = parameters(2,1:4,2);
p_ub = parameters(3,1:4,2);
end

% Make auto fit
options = optimset('Display', 'off');
p_fit_temp = lsqnonlin(residual,p_i,p_lb,p_ub,options);
p_fit(n,2:size(p_fit_temp,2)+1,m) = p_fit_temp;    % plug in current fit values to p_fit matrix

% Plot fits, if specified
if run_op.q_plots == 1
    
gauss_1 = p_fit_temp(1)*exp(-k_eq.^2*p_fit_temp(2)^2);    % Diffusion Gaussian fit
noise = p_fit_temp(3).^2*ss./abs(1-p_fit_temp(4)*exp(-1i*k_eq*ss)).^2;   % noise fit

% Plot fit results
if isnan(iso.P(:,1)) > 0  % if iso = d17o   
    fig(n+length(iso.P(1,:)),'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
elseif iso.P(4,1) > 1    % if iso = dD, plot in different figures
    fig(n+2*length(iso.P(1,:)),'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
else                    % if iso = d18o
    fig(n,'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
end

%semilogy(f_eq,P_eq,'b','LineWidth',1.2)     % plot log-binned spectrum
semilogy(f,G_data,'b','LineWidth',1.2)     % plot original spectrum
hold on
if m ==1    % plot solid lines for actual data and fits
semilogy(f_eq,G(p_fit_temp),'r','LineWidth',1.2)
semilogy(f_eq,gauss_1,'color', [198 120 15] ./ 255,'LineWidth',1.2)
semilogy(f_eq,noise,'k','LineWidth',1.2)
else        % plot dashed lines for uncertainty levels
    semilogy(f_eq,G(p_fit_temp),'r','LineWidth',0.7,'linestyle','--')
    hold on
    semilogy(f_eq,gauss_1,'color', [198 120 15] ./ 255,'LineWidth',0.7,'linestyle','--')
    semilogy(f_eq,noise,'k','LineWidth',0.7,'linestyle','--')
end
xlabel('Frequency [cycles/m]')
ylabel(['PSD [',char(8240),'^2 \cdot m]'])
title(['Age = ' num2str(age(n)) ' yr BP'])
legend('Data','Total Fit', 'Gaussian', 'Noise')
if iso.P(4,1) > 1       % determine if we're scaling axes for dD
    axis([0 100 5*10^-6 5*10^1])
else                    % of if we're scaling axes of d18o or d17o
    axis([0 100 7*10^-7 10^0])
end
set(gca,'YminorTick','off')

end
end


% Double Gaussian Fit
if run_op.q_2gauss == 1

% Define model for fit
G = @(p) p(1)*exp(-k_eq.^2*p(2)^2)...      % Gaussian diffusion function
    + p(3).^2*ss./abs(1-p(4)*exp(-1i*k_eq*ss)).^2 ... % Red noise function
    + p(5)*exp(-k_eq.^2*p(6)^2);               % second gaussian function'
    

% calculate residual of model - data
residual = @(p) log10(G(p)) - log10(P_eq);


%Initial parameters and limits for 2 gaussians
% for holocene spectra
if age(n) <= glacial_start
p_i = parameters(1,1:end,1);     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = parameters(2,1:end,1);
p_ub = parameters(3,1:end,1);
% for glacial spectra
elseif age(n) > glacial_start
p_i = parameters(1,1:end,2);     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = parameters(2,1:end,2);
p_ub = parameters(3,1:end,2);
end

% Make auto fit
options = optimset('Display', 'off');
p_fit_temp = lsqnonlin(residual,p_i,p_lb,p_ub,options);
p_fit(n,2:size(p_fit_temp,2)+1,m) = p_fit_temp;    % plug in current fit values to p_fit matrix

% Plot fit results, if specified
if run_op.q_plots == 1

gauss_1 = p_fit_temp(1)*exp(-k_eq.^2*p_fit_temp(2)^2);    % Diffusion Gaussian fit
noise = p_fit_temp(3).^2*ss./abs(1-p_fit_temp(4)*exp(-1i*k_eq*ss)).^2;   % noise fit
gauss_2 = p_fit_temp(5)*exp(-k_eq.^2*p_fit_temp(6)^2); % second Gaussian fit

% Make plot
if isnan(iso.P(:,1)) > 0  % if iso = d17o   
    fig(n+length(iso.P(1,:)),'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
elseif iso.P(4,1) > 1    % if iso = dD, plot in different figures
    fig(n+2*length(iso.P(1,:)),'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
else                    % if iso = d18o
    fig(n,'units','inches','width',8.5,'height',5,'font','Helvetica','fontsize',15);
end
% loglog(f_eq,P_eq,'b','LineWidth',1.2)          % plot binned spectrum
loglog(f,G_data,'b','LineWidth',1.2)         % plot un-binned spectrum
hold on

% plot main fits
if m ==1    
loglog(f_eq,G(p_fit_temp),'r','LineWidth',1.2)
hold on
loglog(f_eq,gauss_1,'color', [198 120 15] ./ 255,'LineWidth',1.2)
loglog(f_eq,noise,'k','LineWidth',1.2)
loglog(f_eq,gauss_2,'color', [16 149 70] ./ 255,'LineWidth',1.2)

% plot uncertainty bounds
else        
    loglog(f_eq,G(p_fit_temp),'r','LineWidth',0.7,'linestyle','--')
    hold on
    loglog(f_eq,gauss_1,'color', [198 120 15] ./ 255,'LineWidth',0.7,'linestyle','--')
    loglog(f_eq,noise,'k','LineWidth',0.7,'linestyle','--')
    loglog(f_eq,gauss_2,'color', [16 149 70] ./ 255,'LineWidth',0.7,'linestyle','--')
end

xlabel('Frequency [cycles/m]')
ylabel(['PSD ','[',char(8240),'^2 \cdot m]'])
title([num2str(age(n)) ' ' spec_info.win_units])      % label depth/age of middle of window
legend('Data','Total Fit', 'Gaussian', 'Noise','Second Gaussian')
if iso.P(4,1) > 1       % determine if we're scaling axes for dD
    axis([0 100 5*10^-6 5*10^1])
else                    % of if we're scaling axes of d18o or d17o
    axis([0 100 7*10^-7 10^0])
end
set(gca,'YminorTick','off')
end
end

end
end

end