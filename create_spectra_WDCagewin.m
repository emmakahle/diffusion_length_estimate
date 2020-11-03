function [dD_P_matrix_emma, d18o_P_matrix, d17o_P_matrix,start_yr,start_yr_in,end_yr_in] = create_spectra_WDCagewin(method,datafile,win_len_yr,step,ARu,noise_toggle)

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


oldfolder = cd('/Users/emma/Documents/Research/Data');
data = importdata(datafile);    % import data file
cd(oldfolder)
depth = data(69:680786,1);              % assign depth, only include data within age scale bounds
dD = data(69:680786,2);                 % assign dD
d18o = data(69:680786,3);               % assign d18o
if size(data,2)>3
    d17o = data(69:680786,4);           % assign d17o, if it is there
end

% add noise to data if doing 
if noise_toggle == 1
    % add noise to Data
    d18O_with_noise = d18o + 1.7*rand(size(d18o));
    dD_with_noise = dD + 7*rand(size(dD));
end


% import age scale
oldfolder = cd('/Users/emma/Documents/Research/Data');
WAIS_depth_age = importdata('WD2014_agescale.txt');
cd(oldfolder)

% interp to age scale
diff_age = interp1(WAIS_depth_age.data(:,1),WAIS_depth_age.data(:,2),depth);

% start_yr = 1:step:55000-win_len_yr;     % create start year for each window
% end_yr = start_yr+win_len_yr;                   % create end year for each window
% for n = 1:length(start_yr)      
% start_yr_in(n) = find(diff_age>=start_yr(n),1,'first'); % find index for start of windows (same for depth and age)
% end_yr_in(n) = find(diff_age>=end_yr(n),1,'first');     % find index of end of windows (same for depth and age)
% end

% Create same windows as tyler
start_yr = 440-250:step:22940-250;     % create start year for each window
end_yr = start_yr+win_len_yr;                   % create end year for each window
for n = 1:length(start_yr)      
start_yr_in(n) = find(diff_age>=start_yr(n),1,'first'); % find index for start of windows (same for depth and age)
end_yr_in(n) = find(diff_age>=end_yr(n),1,'first');     % find index of end of windows (same for depth and age)
end

middle_yr = sort(repmat(start_yr+(win_len_yr/2-1),1,4));              % find middle meter of each window, repeating 4 times

% % use window from 450-500m depth for paper
% start_m_in = find(depth >= 450,1,'first');
% end_m_in = find(depth >= 500,1,'first');


for n = 1:length(start_yr_in)
    x = dD(start_yr_in(n): end_yr_in(n)); % trims window
    if noise_toggle == 1
        x = dD_with_noise(start_yr_in(n): end_yr_in(n)); % trims window
    end
% for n = 1
%     x = dD(start_m_in:end_m_in); % window for paper 450-500m
    if method == 'MEM'
    %[pxx,fx,pxxc] = pburg(x,ARu,[],1/.005,'ConfidenceLevel',.95);   % calc MEM power, frequency, and 95% confidence
    [pxx,fx] = MEM_1974(x,100,0.005);
%     
%     % add noise to data if doing 
%     if noise_toggle == 1
%     % add noise to Data
%     pxx = pxx + MEM_1974(7*rand(size(dD)),100,1);
%     end
    
    MEM_dD_P_matrix(:,4*n-3) = fx;    % assign results
    MEM_dD_P_matrix(:,4*n-2) = pxx*100;     % make comparable to MTM
    MEM_dD_P_matrix(:,4*n-1) = pxxc(:,1);
    MEM_dD_P_matrix(:,4*n) = pxxc(:,2);
    end
    if method == 'MTM'
    [P,F,C] = pmtmPH(x,.005,10,0,2000);
    
%       % add noise to data if doing 
%     if noise_toggle == 1
%     % add noise to Data
%     P = P + pmtmPH(7*rand(size(dD)),1,10,0,2000);
%     end
    
    %P = P*.005;     % adjust for possible units error in PSD (multiply by same size (m))
    cu=P.*C(:,1);    %get spectrum estimate from confidence level
    cl=P.*C(:,2);
    f_start = find(F==2);
    F=F(f_start:end);
    P=P(f_start:end);
    cu=cu(f_start:end);
    cl=cl(f_start:end);
    MTM_dD_P_matrix_emma(:,4*n-3) = F;    % assign results
    MTM_dD_P_matrix_emma(:,4*n-2) = P;
    MTM_dD_P_matrix_emma(:,4*n-1) = cu;
    MTM_dD_P_matrix_emma(:,4*n) = cl;
    end 
end

if method == 'MEM'
    dD_P_matrix_emma = [middle_yr; MEM_dD_P_matrix];      % Put depths (middle of window) at top of matrix
end
if method == 'MTM'
    dD_P_matrix_emma = [middle_yr; MTM_dD_P_matrix_emma];      % Put depths (middle of window) at top of matrix
end

% % for paper
% if method == 'MEM'
%     dD_P_matrix_emma = [repmat(425,1,4); MEM_dD_P_matrix];      % Put depths (middle of window) at top of matrix
% end
% if method == 'MTM'
%     dD_P_matrix_emma = [repmat(425,1,4); MTM_dD_P_matrix_emma];      % Put depths (middle of window) at top of matrix
% end


% Calculate d18o spectra
for n = 1:length(start_yr_in)
    x = d18o(start_yr_in(n): end_yr_in(n)); % trims window
    if noise_toggle == 1
        x = d18O_with_noise(start_yr_in(n): end_yr_in(n)); % trims window
    end
    if method == 'MEM'
    [pxx,fx,pxxc] = pburg(x,ARu,[],1/.005,'ConfidenceLevel',.95);   % calc MEM power, frequency, and 95% confidence
        
%     % add noise to data if doing 
%     if noise_toggle == 1
%     % add noise to Data
%     pxx = pxx + MEM_1974(1.7*rand(size(d18o)),100,1);
%     end
    
    MEM_d18o_P_matrix(:,4*n-3) = fx;    % assign results
    MEM_d18o_P_matrix(:,4*n-2) = pxx;
    MEM_d18o_P_matrix(:,4*n-1) = pxxc(:,1);
    MEM_d18o_P_matrix(:,4*n) = pxxc(:,2);
    end
    if method == 'MTM'
    [P,F,C] = pmtmPH(x,.005,10,0,2000);
    
%     % add noise to data if doing 
%     if noise_toggle == 1
%     % add noise to Data
%     P = P + pmtmPH(1.7*rand(size(d18o)),1,10,0,2000);
%     end
%    
%     
    P = P*.005;     % adjust for possible units error in PSD (multiply by same size (m))
    cu=P.*C(:,1);    %get spectrum estimate from confidence level
    cl=P.*C(:,2);
    MTM_d18o_P_matrix(:,4*n-3) = F(6:end);    % assign results
    MTM_d18o_P_matrix(:,4*n-2) = P(6:end);
    MTM_d18o_P_matrix(:,4*n-1) = cu(6:end);
    MTM_d18o_P_matrix(:,4*n) = cl(6:end);
    end     
end

if method == 'MEM'
    d18o_P_matrix = [middle_yr; MEM_d18o_P_matrix];      % Put depths (middle of window) at top of matrix
end
if method == 'MTM'
    d18o_P_matrix = [middle_yr; MTM_d18o_P_matrix];      % Put depths (middle of window) at top of matrix
end

if  size(data,2)>3
% Calculate d17o spectra
for n = 1:length(start_yr_in)
    x = d17o(start_yr_in(n): end_yr_in(n)); % trims window
    if method == 'MEM'
    [pxx,fx,pxxc] = pburg(x,ARu,[],1/.005,'ConfidenceLevel',.95);   % calc MEM power, frequency, and 95% confidence
    MEM_d17o_P_matrix(:,4*n-3) = fx;    % assign results
    MEM_d17o_P_matrix(:,4*n-2) = pxx;
    MEM_d17o_P_matrix(:,4*n-1) = pxxc(:,1);
    MEM_d17o_P_matrix(:,4*n) = pxxc(:,2);
    end
    if method == 'MTM'
    [P,F,C] = pmtmPH(x,.005,10,0,2000);
    P = P*.005;     % adjust for possible units error in PSD (multiply by same size (m))
    cu=P.*C(:,1);    %get spectrum estimate from confidence level
    cl=P.*C(:,2);
    MTM_d17o_P_matrix(:,4*n-3) = F(6:end);    % assign results
    MTM_d17o_P_matrix(:,4*n-2) = P(6:end);
    MTM_d17o_P_matrix(:,4*n-1) = cu(6:end);
    MTM_d17o_P_matrix(:,4*n) = cl(6:end);
    end     
end


if method == 'MEM'
    d17o_P_matrix = [middle_yr; MEM_d17o_P_matrix];      % Put depths (middle of window) at top of matrix
end
if method == 'MTM'
    d17o_P_matrix = [middle_yr; MTM_d17o_P_matrix];      % Put depths (middle of window) at top of matrix
end
end

% % Save matrices of spectral data
% dlmwrite(['SP_spectra_dD_' num2str(win_len_m) 'm_' num2str(step_m) 'm.txt'], dD_P_matrix, 'delimiter', '\t', 'precision', 7)
% dlmwrite(['SP_spectra_d18o_' num2str(win_len_m) 'm_' num2str(step_m) 'm.txt'], d18o_P_matrix, 'delimiter', '\t', 'precision', 7)
% %dlmwrite(['MTM_SP_spectra_dD_' num2str(win_len_m) 'm_' num2str(step_m) 'm.txt'], MTM_dD_P_matrix, 'delimiter', '\t', 'precision', 7)
% 
% savename_dD = ['SP_spectra_dD_' num2str(win_len_m) 'm_' num2str(step_m) 'm.txt'];
% savename_d18o = ['SP_spectra_d18o_' num2str(win_len_m) 'm_' num2str(step_m) 'm.txt'];
end