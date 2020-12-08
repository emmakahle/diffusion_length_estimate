%% SPC diffusion project code

% Author: Emma C. Kahle
% Last Update: 11/3/2020
% Details in corresponding publication: 
% Kahle, E. C., Holme, C., Jones, T. R., Gkinis, V., & Steig, E. J. (2018). 
% A Generalized Approach to Estimating Diffusion Length of Stable Water 
% Isotopes From Ice-Core Data. Journal of Geophysical Research: 
% Earth Surface, 123(10), 2377– 2391. https://doi.org/10.1029/2018JF004764

% This is the main code to run to create diffusion length estimates from
% ice core water isotope data. This can be utilized with other water
% isotope data sets as well.

%% Set up run options

% run options (1 turns a setting on, 0 turns a setting off)
run_op.q_1gauss = 0;   % use single gaussian parameterization
run_op.q_2gauss = 1;   % use double gaussian parameterization
run_op.q_plots = 0;    % make plots of each spectrum fit (1 cues plots, 0 does not plot)
run_op.remove_peak = 0;    % remove annual peak from spectra
run_op.uncertainty = 1;    % calculate uncertainty on sigma estimate
run_op.noise_add = 0;      % also calculate noise-adding technique

% save options
save_spectra = 0; % save spectra info for each isotope
save_sigma = 0;     % save diffusion lengths for each isotope

% spectra options
ARu = 100;      % number of poles for MEM analysis
method = 'MEM'; % choose 'MEM' or 'MTM'
%spec_info.win_units = 'yr';    % choose to make windows either in constant time 'yr' or constant depth 'm'
%spec_info.win_len = 250;    % in whichever unit you specify above (either years or meters)
%spec_info.step_len = 250;       % in whichever unit you specify above (either years or meters)
spec_info.win_units = 'm';    % choose to make windows either in constant time 'yr' or constant depth 'm'
spec_info.win_len = 10;    % in whichever unit you specify above (either years or meters)
spec_info.step_len = 10;       % in whichever unit you specify above (either years or meters)


% Define Holocene and Glacial regimes - can use different schemes for
% different time periods
glacial_start = 10000; % year glacial regime begins

%% Import/create spectra
% Identify data files to with water isotope data. Data must be evenly
% sampled at half-cm intervals.

% Set parameters for creating spectra
spec_info.save = save_spectra;
noise_toggle = 0;

% names of datafiles
datafile17 = 'Version1_201520162017_5mm_interp.mat';
datafile = 'SPC_ISO_halfcm_INSTAAR_linearNaN_depth_d18o_dD.txt';

% load water isotope data
data17 = load(datafile17);    % import data file
data_i = importdata(datafile);    % import data file

% interp d17O data to rest of data
d17O = interp1(data17.depthi, data17.d17Oi, data_i(:,1));

% pull all data in one matrix
data = [data_i(:,1),data_i(:,2),data_i(:,3),d17O];

% Create spectra with given parameters
[iso_spec] = create_spectra(method,data,spec_info,ARu,noise_toggle);


%% Autofit 

% Set fit parameters (for d18o)
% HOLOCENE
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.1 0 .01 0 0 0];
p_ub = [100 .12 10 .5 7e-04 .015];  
parameters_d18o_hol = [p_i;p_lb;p_ub];
% GLACIAL
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.1 0 .01 0 0 0];
p_ub = [100 .12 10 .5 5e-04 .015];  
parameters_d18o_glac = [p_i;p_lb;p_ub];

parameters_d18o(:,:,1) = parameters_d18o_hol;
parameters_d18o(:,:,2) = parameters_d18o_glac;

[p_fit_d18o] = auto_fit(iso_spec.d18o, parameters_d18o, run_op, spec_info, glacial_start);


% Set fit parameters (for d17o)
% HOLOCENE
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.01 0 .01 0 9e-05 0];
p_ub = [100 .12 10 .1 7e-04 .015];  
parameters_d17o_hol = [p_i;p_lb;p_ub];
% GLACIAL
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.01 0 .01 0 9e-05 0];
p_ub = [100 .12 10 .1 2e-04 .015];  
parameters_d17o_glac = [p_i;p_lb;p_ub];

parameters_d17o(:,:,1) = parameters_d17o_hol;
parameters_d17o(:,:,2) = parameters_d17o_glac;

[p_fit_d17o] = auto_fit(iso_spec.d17o, parameters_d17o, run_op, spec_info, glacial_start);

% Set fit parameters (for dD) 
% HOLOCENE
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.1 0 .01 0 1e-4 0];
p_ub = [100 .12 10 .5 1e-01 .015];  
parameters_dD_hol = [p_i;p_lb;p_ub];
% GLACIAL
p_i = [1 .06 2.5e-05 .03 3e-4 .007];     % [P_1, sigma_1, sigma_noise, ar1, P_2, sigma_2]
p_lb = [.1 0 .01 0 1e-4 0];
p_ub = [100 .12 10 .5 1e-02 .015];  
parameters_dD_glac = [p_i;p_lb;p_ub];

parameters_dD(:,:,1) = parameters_dD_hol;
parameters_dD(:,:,2) = parameters_dD_glac;

% run fit
[p_fit_dD] = auto_fit(iso_spec.dD, parameters_dD, run_op, spec_info, glacial_start);


%% Run noise_adding technique if specified at top
if run_op.noise_add == 1
    noise_toggle = 1;
    
    % create spectra with noise added
    [dD_NA, d18o_NA,~,win_params] = create_spectra_WDCagewin(method,datafile,win_len,step_len,ARu,noise_toggle);
    
    % autofit
    % d18o
    p_i = [1 .06 2.5e-05 .03];     % [P_1, sigma_1, sigma_noise, ar1]
    p_lb = [.1 0 .01 0];
    p_ub = [100 .12 10 .5];  
    parameters_d18o_NA = [p_i;p_lb;p_ub];
    [p_fit_d18o_NA] = auto_fit(d18o_NA, parameters_d18o_NA, 1, 0, q_plots, remove_peak, uncertainty);

    % dD
    p_i = [10 .06 5 .1];     % [P_1, sigma_1, sigma_noise, ar1]
    p_lb = [.1 0 .01 0];
    p_ub = [10^4 .12 100 .5];  
    parameters_dD_NA = [p_i;p_lb;p_ub];
    [p_fit_dD_NA] = auto_fit_compare_Tyler(dD_NA, parameters_dD_NA, 1, 0, q_plots, remove_peak, uncertainty);
end


%% Plot results

% dD fits
auto_age_dD = iso_spec.dD.P(1,:);
auto_sigmaD = p_fit_dD(:,3,1);
auto_sigma_upD = p_fit_dD(:,3,2);
auto_sigma_dnD = p_fit_dD(:,3,3);

% auto_age_NA = dD_NA.P(1,:);
% auto_sigmaD_NA = p_fit_dD_NA(:,3,1);
% auto_sigma_upD_NA = p_fit_dD_NA(:,3,2);
% auto_sigma_dnD_NA = p_fit_dD_NA(:,3,3);

% d18o fits
auto_age_d18o = iso_spec.d18o.P(1,:);
auto_sigma18o = p_fit_d18o(:,3,1);
auto_sigma_up18o = p_fit_d18o(:,3,2);
auto_sigma_dn18o = p_fit_d18o(:,3,3);

% auto_age_NA = d18o_NA.P(1,:);
% auto_sigma18o_NA = p_fit_d18o_NA(:,3,1);
% auto_sigma_up18o_NA = p_fit_d18o_NA(:,3,2);
% auto_sigma_dn18o_NA = p_fit_d18o_NA(:,3,3);

% d17o fits
d17o_start = find(p_fit_d17o(:,3,1)>0, 1, 'first'); % find start of d17o data
auto_age_d17o = iso_spec.d17o.P(1,d17o_start:end);
auto_sigma17o = p_fit_d17o(d17o_start:end,3,1);
auto_sigma_up17o = p_fit_d17o(d17o_start:end,3,2);
auto_sigma_dn17o = p_fit_d17o(d17o_start:end,3,3);

% adjust uncertainty ranges to account for any upper bounds that actually
% came out as lower
for ii = 1:size(auto_sigma17o,1)
    if auto_sigma_up17o(ii)>auto_sigma17o(ii)   % if upper bound is actually lower
        auto_sigma_up17o(ii) = auto_sigma17o(ii) - (auto_sigma17o(ii-1)-auto_sigma_up17o(ii-1));   % set it to uncertainty of previous window
    end
end

% auto_age_NA = d17o_NA.P(1,:);
% auto_sigma17o_NA = p_fit_d17o_NA(:,3,1);
% auto_sigma_up17o_NA = p_fit_d17o_NA(:,3,2);
% auto_sigma_dn17o_NA = p_fit_d17o_NA(:,3,3);


% Plot d18o and dD for NA and DG - 1 panel figure
fig('units','inches','width',10,'height',6,'font','Helvetica','fontsize',16);
if run_op.uncertainty == 1
ax1 = shadedErrorBar(auto_age_d18o/1000,auto_sigma18o,[auto_sigma18o-auto_sigma_dn18o auto_sigma_up18o-auto_sigma18o]','-b',0.5);   
% hold on
% ax2 = shadedErrorBar(auto_age_NA/1000,auto_sigma18o_NA,[auto_sigma18o_NA-auto_sigma_dn18o_NA auto_sigma_up18o_NA-auto_sigma18o_NA]','-ok',0.5);   
elseif run_op.uncertainty == 0
ax1 = plot(auto_age_d18o/1000,auto_sigma18o,'b','LineWidth',1);
% ax2 = plot(auto_age_NA/1000,auto_sigma18o_NA);
end
% if run_op.uncertainty ==1
% elseif run_op.uncertainty ==0
% end

% d17o
hold on
if run_op.uncertainty == 1
ax2 = shadedErrorBar(auto_age_d17o/1000,auto_sigma17o,[auto_sigma17o-auto_sigma_dn17o auto_sigma_up17o-auto_sigma17o]','-c',0.5);   
% hold on
% ax2 = shadedErrorBar(auto_age_NA/1000,auto_sigma17o_NA,[auto_sigmad17o_NA-auto_sigma_dn17o_NA auto_sigma_up17o_NA-auto_sigma17o_NA]','-ok',0.5);   
elseif run_op.uncertainty == 0
ax2 = plot(auto_age_d17o/1000,auto_sigma17o,'-c','LineWidth',1);
% ax2 = plot(auto_age_NA/1000,auto_sigma17o_NA);
end


% dD
if run_op.uncertainty == 1
ax3 = shadedErrorBar(auto_age_dD/1000,auto_sigmaD,[auto_sigmaD-auto_sigma_dnD auto_sigma_upD-auto_sigmaD]','-r',0.5);   
% hold on
% ax2 = shadedErrorBar(auto_age_NA/1000,auto_sigmaD_NA,[auto_sigmaD_NA-auto_sigma_dnD_NA auto_sigma_upD_NA-auto_sigmaD_NA]','-ok',0.5);   
elseif run_op.uncertainty == 0
ax3 = plot(auto_age_dD/1000,auto_sigmaD,'-r','LineWidth',1);
% ax2 = plot(auto_age_NA/1000,auto_sigmaD_NA);
end
xlabel('Age [ka]')
ylabel('Diffusion Length [m]')

legend('\delta^{18}O','\delta^{17}O','\delta D')

