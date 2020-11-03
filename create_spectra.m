function [iso_spec,win_params] = create_spectra(method,data,spec_info,ARu,noise_toggle)

%% Create spectra

% Create spectra of isotope data with given parameters

% Inputs
% Method: MEM or MTM
% Datafile: file name from which to import ice core data. Column one is
% depth, column 2 is dD, column 3 is d18o.
% win_len: length of desired window for spectra in meters
% step: length of step in between windows in meters
% ARu: Order of autoregressive model for MEM estimate

% Outputs
% dD_P_matrix: Matrix oriented column-wise, containing frequency, power,
% and upper and lower confidence bounds at 95% for each calculated window.
% Width of matrix should be 4 times the number of windows calculated.
% d18o_P_matrix: same as above for d18o



% depth = data.depthi;              % assign depth, only include data within age scale bounds
% dD = data.dDi;                 % assign dD
% d18o = data.d18Oi;               % assign d18o
% d17o = data.d17Oi;
depth = data(:,1);
dD = data(:,3);
d18o = data(:,2);
d17o = data(:,4);

% add noise to data if doing 
if noise_toggle == 1
    % add noise to Data
    d18O_with_noise = d18o + 1.7*rand(size(d18o));
    dD_with_noise = dD + 7*rand(size(dD));
end


% import age scale
oldfolder = cd('/Users/emma/Documents/Research/SPICE');
SPC_depth_age = importdata('SP14-02_agescale_depth_ice_gas.txt');
cd(oldfolder)

if spec_info.win_units == 'yr'
    % interp to age scale
    age = interp1(SPC_depth_age(:,1),SPC_depth_age(:,2),depth);
    age((isnan(age))) = [];     % remove nan values

    % create windows
    start_yr = age(1):spec_info.step_len:age(end)-spec_info.win_len;     % create start year for each window
    end_yr = start_yr+spec_info.win_len;                   % create end year for each window
    middle = sort(repmat(start_yr+(spec_info.win_len/2-1),1,4));              % find middle meter of each window, repeating 4 times
    for n = 1:length(start_yr)      
    start_in(n) = find(age>=start_yr(n),1,'first'); % find index for start of windows (same for depth and age)
    end_in(n) = find(age>=end_yr(n),1,'first');     % find index of end of windows (same for depth and age)
    end
elseif spec_info.win_units == 'm'
    start_m = depth(1):spec_info.step_len:depth(end)-spec_info.win_len;
    end_m = start_m+spec_info.win_len;
    middle = sort(repmat(start_m+(spec_info.win_len/2-1),1,4));              % find middle meter/year of each window, repeating 4 times
    for n = 1:length(start_m)      
    start_in(n) = find(depth>=start_m(n),1,'first'); % find index for start of windows (same for depth and age)
    end_in(n) = find(depth>=end_m(n),1,'first');     % find index of end of windows (same for depth and age)
    end
else
    print('Error: select either y or m for win_units')
end

%%%%%%%%%% can i get rid of this?
% assign data structure to hold window parameters
% win_params.start_yr = start_yr;
% win_params.start_yr_in = start_in;
% win_params.end_yr_in = end_in;

%% Calculate spectra
% input into large matrices for each isotope
for n = 1:length(start_in)
    x = dD(start_in(n): end_in(n)); % trims window
    if noise_toggle == 1
        x = dD_with_noise(start_in(n): end_in(n)); % trims window
    end
    if method == 'MEM'
    [pxx,fx,pxxc] = pburg(x,ARu,[],1/.005,'ConfidenceLevel',.68);   % calc MEM power, frequency, and 95% confidence
%     tf = isempty(pxxc); % flag when no confidence level is calculated
    %[pxx,fx] = MEM_1974(x,100,0.005);  % christian's code does not output
    %confidence levels
    full_dD_P_matrix(:,4*n-3) = fx;    % assign results
    full_dD_P_matrix(:,4*n-2) = pxx;
%     if tf == 0  % if confidence level is calculated, fill with pxxc
    full_dD_P_matrix(:,4*n-1) = pxxc(:,1);
    full_dD_P_matrix(:,4*n) = pxxc(:,2);
%     else    % if no confidence level is calculated, fill with nan
%     full_dD_P_matrix(:,4*n-1) = NaN(size(fx));
%     full_dD_P_matrix(:,4*n) = NaN(size(fx));
%     end
    end
    if method == 'MTM'
    [P,F,C] = pmtmPH(x,.005,10,0,2000);
    cu=P.*C(:,1);    %get spectrum estimate from confidence level
    cl=P.*C(:,2);
    f_start = find(F==2);
    full_dD_P_matrix(:,4*n-3) = F(f_start:end);    % assign results
    full_dD_P_matrix(:,4*n-2) = P(f_start:end);
    full_dD_P_matrix(:,4*n-1) = cu(f_start:end);
    full_dD_P_matrix(:,4*n) = cl(f_start:end);
    end 
end

dD_P_matrix = [middle; full_dD_P_matrix];      % Put depths (middle of window) at top of matrix


% Calculate d18o spectra
for n = 1:length(start_in)
    x = d18o(start_in(n): end_in(n)); % trims window
    if noise_toggle == 1
        x = d18O_with_noise(start_in(n): end_in(n)); % trims window
    end
    if method == 'MEM'
    [pxx,fx,pxxc] = pburg(x,ARu,[],1/.005,'ConfidenceLevel',.68);   % calc MEM power, frequency, and 95% confidence
    full_d18o_P_matrix(:,4*n-3) = fx;    % assign results
    full_d18o_P_matrix(:,4*n-2) = pxx;
    full_d18o_P_matrix(:,4*n-1) = pxxc(:,1);
    full_d18o_P_matrix(:,4*n) = pxxc(:,2);
    end
    if method == 'MTM'
    [P,F,C] = pmtmPH(x,.005,10,0,2000);
    cu=P.*C(:,1);    %get spectrum estimate from confidence level
    cl=P.*C(:,2);
    full_d18o_P_matrix(:,4*n-3) = F(6:end);    % assign results
    full_d18o_P_matrix(:,4*n-2) = P(6:end);
    full_d18o_P_matrix(:,4*n-1) = cu(6:end);
    full_d18o_P_matrix(:,4*n) = cl(6:end);
    end     
end

d18o_P_matrix = [middle; full_d18o_P_matrix];      % Put depths (middle of window) at top of matrix

% Calculate d17o spectra
for n = 1:length(start_in)
    x = d17o(start_in(n): end_in(n)); % trims window
    if sum(isnan(x))==0
        if method == 'MEM'
        [pxx,fx,pxxc] = pburg(x,ARu,[],1/.005,'ConfidenceLevel',.68);   % calc MEM power, frequency, and 95% confidence
        full_d17o_P_matrix(:,4*n-3) = fx;    % assign results
        full_d17o_P_matrix(:,4*n-2) = pxx;
        full_d17o_P_matrix(:,4*n-1) = pxxc(:,1);
        full_d17o_P_matrix(:,4*n) = pxxc(:,2);
        end
        if method == 'MTM'
        [P,F,C] = pmtmPH(x,.005,10,0,2000);
        cu=P.*C(:,1);    %get spectrum estimate from confidence level
        cl=P.*C(:,2);
        full_d17o_P_matrix(:,4*n-3) = F(6:end);    % assign results
        full_d17o_P_matrix(:,4*n-2) = P(6:end);
        full_d17o_P_matrix(:,4*n-1) = cu(6:end);
        full_d17o_P_matrix(:,4*n) = cl(6:end);
        end     
    else
        full_d17o_P_matrix = nan(size(full_d18o_P_matrix));    % assign nan if d17o data doesn't exist
    end
end

d17o_P_matrix = [middle; full_d17o_P_matrix];      % Put depths (middle of window) at top of matrix


%% Assign data into structures for use
% inputs matrices with spectral estimates for dD and d18o, and outputs
% structures for dD and d18o that contain frequency, power, and upper and
% lower 95% confidence level estimates for power

% Assign dD data
dD_spec.f = dD_P_matrix(:,1); % assign frequency vector for all power vectors
dD_spec.P = dD_P_matrix(:, 2:4:size(dD_P_matrix,2));    % assign power vectors in one matrix
dD_spec.cu = dD_P_matrix(:, 3:4:size(dD_P_matrix,2));   % assign 95% confidence estimate
dD_spec.cl = dD_P_matrix(:, 4:4:size(dD_P_matrix,2));

% Assign d18o data
d18o_spec.f = d18o_P_matrix(:,1); % assign frequency vector for all power vectors
d18o_spec.P = d18o_P_matrix(:, 2:4:size(d18o_P_matrix,2));    % assign power vectors in one matrix
d18o_spec.cu = d18o_P_matrix(:, 3:4:size(d18o_P_matrix,2));   % assign 95% confidence estimate
d18o_spec.cl = d18o_P_matrix(:, 4:4:size(d18o_P_matrix,2));

% Assign d17o data
d17o_spec.f = d17o_P_matrix(:,end-3); % assign frequency vector for all power vectors
d17o_spec.P = d17o_P_matrix(:, 2:4:size(d17o_P_matrix,2));    % assign power vectors in one matrix
d17o_spec.cu = d17o_P_matrix(:, 3:4:size(d17o_P_matrix,2));   % assign 95% confidence estimate
d17o_spec.cl = d17o_P_matrix(:, 4:4:size(d17o_P_matrix,2));

% Create structure for all isotope spectra information
iso_spec.dD = dD_spec;
iso_spec.d18o = d18o_spec;
iso_spec.d17o = d17o_spec;


%% Save matrices of spectral data (if specified)
if spec_info.save == 1
    dlmwrite(['SP_spectra_dD_' num2str(spec_info.win_len) spec_info.win_units num2str(spec_info.step_len) spec_info.win_units '.txt'], dD_P_matrix, 'delimiter', '\t', 'precision', 7)
    dlmwrite(['SP_spectra_d18o_' num2str(spec_info.win_len) spec_info.win_units num2str(spec_info.step_len) spec_info.win_units '.txt'], d18o_P_matrix, 'delimiter', '\t', 'precision', 7)
    dlmwrite(['SP_spectra_d17o_' num2str(spec_info.win_len) spec_info.win_units num2str(spec_info.step_len) spec_info.win_units '.txt'], d17o_P_matrix, 'delimiter', '\t', 'precision', 7)
end

end