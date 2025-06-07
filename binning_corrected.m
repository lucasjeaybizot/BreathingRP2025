% ------------------------------------------------------------------------
% This script loads and processes EEG and respiration phase data from three
% sources: original, new, and simulated. It bins EEG amplitude according to
% our corrected binning method, and generates permutation distributions.
% The corrected method groups trials according to the phase the button
% press occured in Modulation Index permutation statistices are also
% computed.
% The circ_stat toolbox is required.
% ------------------------------------------------------------------------

%% Fix random seed for reproducibility
rng(2024)

%% --------------------------- Original Data --------------------------- %%
clear
files_original                      = dir('data\original\processed\data_EEG_Suj*_E.mat');
originaldata_correctedbinning       = struct();

for participant = 1:length(files_original)
    % Load data
    load([files_original(participant).folder '\' sprintf('data_EEG_Suj%02d_E.mat', participant)])
    load([files_original(participant).folder '\' sprintf('Resp_phase_Suj%02d_E.mat', participant)])
    
    % Eval EEG and Resp
    dat     = eval(sprintf('data_EEG_Suj%02d_E', participant));
    resp    = eval(sprintf('Resp_phase_Suj%02d_E', participant));

    % Clear memory
    clear(sprintf('data_EEG_Suj%02d_E', participant), sprintf('Resp_phase_Suj%02d_E', participant));
    
    % Keep 2 s before press for RP amplitude
    dat.trial   = cellfun(@(x) x(1024:2048),dat.trial,'UniformOutput',false);
    RP_amp      = mean(reshape(cell2mat(dat.trial'), [size(resp,1) length(1024:2048)]),2);

    % Bin according to Resp phase at time of button press
    Resp_phase_press = resp(:,2048);
    
    bin_ID              = 0;
    mean_amp            = nan(1,6);

    for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
        bin_ID = bin_ID + 1;
        I = find(Resp_phase_press > phase_bin & Resp_phase_press <= phase_bin + pi/3);
        mean_amp(bin_ID) = mean(RP_amp(I),'all');
    end

    % Get permutations
    perm_amp = nan(1000,6);
    
    for i = 1:1000
        perm_idx = randperm(length(Resp_phase_press));
        bin_ID = 0;

        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find(Resp_phase_press(perm_idx) > phase_bin & Resp_phase_press(perm_idx) <= phase_bin + pi/3);
            perm_amp(i,bin_ID) = mean(RP_amp(I),'all');
        end
    end
    
    % Save files
    originaldata_correctedbinning(participant).Resp_phase_press   = Resp_phase_press;        % Used for rose plot
    originaldata_correctedbinning(participant).mean_amp           = mean_amp;                % Used for Bayesian test and MI test
    originaldata_correctedbinning(participant).perm_amp           = perm_amp;                % Used for MI test
    originaldata_correctedbinning(participant).mean_amp_sumnorm   = mean_amp ./ sum(mean_amp); % sum norm from Park et al. 2020 [applied to whole participant vs. per trial due to change of binning]
end

fprintf('Original data binning done. Starting permutation statistics ... \n')

%% Permutation statistics - Original Data

for participant = 1:length(originaldata_correctedbinning)  
    temp_amp        = originaldata_correctedbinning(participant).mean_amp + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
    n               = length(temp_amp);
    p               = temp_amp / sum(temp_amp);
    H               = -sum(p.*log(p));
    originaldata_correctedbinning(participant).MI = (log(n) - H) / log(n);
    
    originaldata_correctedbinning(participant).MI_perm = nan(1, 1000);

    % Permutations
    for i = 1:1000
        temp_amp = originaldata_correctedbinning(participant).perm_amp(i,:) + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
        n = length(temp_amp);
        p = temp_amp / sum(temp_amp);
        H = -sum(p.*log(p));
        originaldata_correctedbinning(participant).MI_perm(i) = (log(n) - H) / log(n);
    end

    % Individual significance (not used)
    originaldata_correctedbinning(participant).pval_MI = sum(originaldata_correctedbinning(participant).MI_perm > originaldata_correctedbinning(participant).MI) / 1000;
end

if any(isnan([originaldata_correctedbinning(:).MI_perm]),'all')
    error('Carefull histogram will not be accurate: N/A values in MI_perm')
end

sum_MI_perm_vals    = sum(reshape([originaldata_correctedbinning.MI_perm], [], length(originaldata_correctedbinning))',1);
sum_MI_vals         = sum(reshape([originaldata_correctedbinning.MI], [], length(originaldata_correctedbinning))',1);
originaldata_correctedbinning(1).pval_MI_total = (sum(sum_MI_perm_vals > sum_MI_vals) + 1) / 1001;

% Outlier rejection 3*SD
perm_data       = reshape([originaldata_correctedbinning.MI_perm], [], length(originaldata_correctedbinning))';
orig_data       = reshape([originaldata_correctedbinning.MI], [], length(originaldata_correctedbinning))';
reject_thres    = 3*std(orig_data);
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
originaldata_correctedbinning(1).pval_MI_total_outlier_3sd = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

% Outlier rejection median-IQR method
outlier_iqr     = prctile(orig_data,75)-prctile(orig_data,25);
reject_thres    = prctile(orig_data,75)+1.5*outlier_iqr;
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
originaldata_correctedbinning(1).pval_MI_total_outlier_iqr = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;


save('results\corrected_binning\originaldata_correctedbinning.mat', 'originaldata_correctedbinning')
fprintf('Original data binning done and saved \n')

%% ----------------------------- New Data ------------------------------ %%
clear
files_new                           = dir('data\new\processed\P*_data_cube.mat');
newdata_correctedbinning            = struct();

for participant = 1:length(files_new)
    % Load data
    load([files_new(participant).folder '\' sprintf('P%02d_data_cube.mat', participant)])

    % Eval EEG and Resp
    dat     = squeeze(data_cube(1,:,:));
    resp    = squeeze(data_cube(2,:,:));

    % Clear memory
    clear('data_cube');
    
    % Keep 2 s before press 
    RP_amp      = mean(dat(1024:2048,:), 1);

    % Bin according to Resp phase at time of button press
    Resp_phase_press = resp(2048,:); 
    
    % .eeg file saves the phase inversed - +pi needed to correct for that (check with respiration_test.eeg file).
    Resp_phase_press_tp = Resp_phase_press;
    Resp_phase_press(Resp_phase_press_tp > 0) = Resp_phase_press_tp(Resp_phase_press_tp > 0) - pi;
    Resp_phase_press(Resp_phase_press_tp < 0) = Resp_phase_press_tp(Resp_phase_press_tp < 0) + pi;

    bin_ID = 0;
    mean_amp = nan(1,6);
    for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3] % [0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3]%
        bin_ID = bin_ID + 1;
        I = find(Resp_phase_press > phase_bin & Resp_phase_press <= phase_bin + pi/3);
        mean_amp(bin_ID) = mean(RP_amp(I),'all');
    end

    % Get permutations
    perm_amp = nan(1000,6);
    
    for i = 1:1000
        perm_idx = randperm(length(Resp_phase_press));
        bin_ID = 0;

        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find(Resp_phase_press(perm_idx) > phase_bin & Resp_phase_press(perm_idx) <= phase_bin + pi/3);
            perm_amp(i,bin_ID) = mean(RP_amp(I),'all');
        end
    end
    
    % Save files
    newdata_correctedbinning(participant).Resp_phase_press   = Resp_phase_press;        % Used for rose plot
    newdata_correctedbinning(participant).mean_amp           = mean_amp;                % Used for Bayesian test and MI test
    newdata_correctedbinning(participant).perm_amp           = perm_amp;                % Used for MI test
    newdata_correctedbinning(participant).mean_amp_sumnorm   = mean_amp ./ sum(mean_amp); % sum norm from Park et al. 2020 [applied to whole participant vs. per trial due to change of binning]
end

fprintf('New data binning done. Starting permutation statistics ... \n')

%% Permutation statistics - New Data

for participant = 1:length(newdata_correctedbinning)  
    temp_amp = newdata_correctedbinning(participant).mean_amp + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
    n = length(temp_amp);
    p = temp_amp / sum(temp_amp);
    H = -sum(p.*log(p));
    newdata_correctedbinning(participant).MI = (log(n) - H) / log(n);
    
    newdata_correctedbinning(participant).MI_perm = nan(1, 1000);

    % Permutations
    for i = 1:1000
        temp_amp = newdata_correctedbinning(participant).perm_amp(i,:) + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
        n = length(temp_amp);
        p = temp_amp / sum(temp_amp);
        H = -sum(p.*log(p));
        newdata_correctedbinning(participant).MI_perm(i) = (log(n) - H) / log(n);
    end

    % Individual significance (% Not used?)
    newdata_correctedbinning(participant).pval_MI = sum(newdata_correctedbinning(participant).MI_perm > newdata_correctedbinning(participant).MI) / 1000; % (not used?)
end

sum_MI_perm_vals    = sum(reshape([newdata_correctedbinning.MI_perm], [], length(newdata_correctedbinning))',1);
sum_MI_vals         = sum(reshape([newdata_correctedbinning.MI], [], length(newdata_correctedbinning))',1);
newdata_correctedbinning(1).pval_MI_total = (sum(sum_MI_perm_vals > sum_MI_vals) + 1) / 1001;

% Outlier rejection
perm_data       = reshape([newdata_correctedbinning.MI_perm], [], length(newdata_correctedbinning))';
orig_data       = reshape([newdata_correctedbinning.MI], [], length(newdata_correctedbinning))';
reject_thres    = 3*std(orig_data);
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
newdata_correctedbinning(1).pval_MI_total_outlier_3sd = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

% Outlier rejection median-IQR method
outlier_iqr     = prctile(orig_data,75)-prctile(orig_data,25);
reject_thres    = prctile(orig_data,75)+1.5*outlier_iqr;
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
newdata_correctedbinning(1).pval_MI_total_outlier_iqr = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

save('results\corrected_binning\newdata_correctedbinning.mat', 'newdata_correctedbinning')
fprintf('New data binning done and saved \n')


%% -------------------------- Simulated Data --------------------------- %%
clear
files_simulated                  = dir('data\simulated\processed\SIMUL*_cube.mat');
simulateddata_correctedbinning    = struct();

for participant = 1:length(files_simulated)
    % Load data
    load([files_simulated(participant).folder '\' sprintf('SIMUL%02d_cube.mat', participant)])

    % Eval EEG and Resp
    dat     = squeeze(data_cube(1,:,:));
    resp    = squeeze(data_cube(2,:,:));

    % Clear memory
    clear('data_cube');
    
    % Keep 2 s before press 
    RP_amp      = mean(dat(1024:2048,:), 1);

    % Bin according to Resp phase at time of button press
    Resp_phase_press = resp(2048,:);
    
    bin_ID = 0;
    mean_amp = nan(1,6);
    for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
        bin_ID = bin_ID + 1;
        I = find(Resp_phase_press >= phase_bin & Resp_phase_press < phase_bin + pi/3);
        mean_amp(bin_ID) = mean(RP_amp(I),'all');
    end

    % Get permutations
    perm_amp = nan(1000,6);
    
    for i = 1:1000
        perm_idx = randperm(length(Resp_phase_press));
        bin_ID = 0;

        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find(Resp_phase_press(perm_idx) >= phase_bin & Resp_phase_press(perm_idx) < phase_bin + pi/3);
            perm_amp(i,bin_ID) = mean(RP_amp(I),'all');
        end
    end
    
    % Save files
    simulateddata_correctedbinning(participant).Resp_phase_press   = Resp_phase_press;        % Used for rose plot
    simulateddata_correctedbinning(participant).mean_amp           = mean_amp;                % Used for Bayesian test and MI test
    simulateddata_correctedbinning(participant).perm_amp           = perm_amp;                % Used for MI test
    simulateddata_correctedbinning(participant).mean_amp_sumnorm   = mean_amp ./ sum(mean_amp); % sum norm from Park et al. 2020 [applied to whole participant vs. per trial due to change of binning]
end

fprintf('Simulated data binning done. Starting permutation statistics ... \n')

%% Permutation statistics - Simulated Data

for participant = 1:length(simulateddata_correctedbinning)  
    temp_amp = simulateddata_correctedbinning(participant).mean_amp + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
    n = length(temp_amp);
    p = temp_amp / sum(temp_amp);
    H = -sum(p.*log(p));
    simulateddata_correctedbinning(participant).MI = (log(n) - H) / log(n);
    
    simulateddata_correctedbinning(participant).MI_perm = nan(1, 1000);

    % Permutations
    for i = 1:1000
        temp_amp = simulateddata_correctedbinning(participant).perm_amp(i,:) + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
        n = length(temp_amp);
        p = temp_amp / sum(temp_amp);
        H = -sum(p.*log(p));
        simulateddata_correctedbinning(participant).MI_perm(i) = (log(n) - H) / log(n);
    end

    % Individual significance
    simulateddata_correctedbinning(participant).pval_MI = sum(simulateddata_correctedbinning(participant).MI_perm > simulateddata_correctedbinning(participant).MI) / 1000; 
end

sum_MI_perm_vals    = sum(reshape([simulateddata_correctedbinning.MI_perm], [], length(simulateddata_correctedbinning))',1);
sum_MI_vals         = sum(reshape([simulateddata_correctedbinning.MI], [], length(simulateddata_correctedbinning))',1);
simulateddata_correctedbinning(1).pval_MI_total = (sum(sum_MI_perm_vals > sum_MI_vals) + 1) / 1001;

% Outlier rejection
perm_data       = reshape([simulateddata_correctedbinning.MI_perm], [], length(simulateddata_correctedbinning))';
orig_data       = reshape([simulateddata_correctedbinning.MI], [], length(simulateddata_correctedbinning))';
reject_thres    = 3*std(orig_data);
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
simulateddata_correctedbinning(1).pval_MI_total_outlier = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

save('results\corrected_binning\simulateddata_correctedbinning.mat', 'simulateddata_correctedbinning')
fprintf('Simulated data binning done and saved \n')