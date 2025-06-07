% ------------------------------------------------------------------------
% This script loads and processes EEG and respiration phase data from three
% sources: original, new, and simulated. It bins EEG amplitude according
% to Park et al. (2020) binning method, and generates permutation
% distributions. Normalized versions of binned amplitudes are also
% computed. Modulation Index permutation statistices are also computed.
% The code herepresent was coded by LJB who had access to Park et al.
% (2020) original code in addition to the description from the published
% manuscript. The circ_stat toolbox is required.
% ------------------------------------------------------------------------

%% Fix random seed for reproducibility
rng(2024)

%% --------------------------- Original Data --------------------------- %%
clear
files_original                  = dir('data\original\processed\data_EEG_Suj*_E.mat');
originaldata_originalbinning    = struct();

for participant = 1:length(files_original)
    % Load data
    load([files_original(participant).folder '\' sprintf('data_EEG_Suj%02d_E.mat', participant)])
    load([files_original(participant).folder '\' sprintf('Resp_phase_Suj%02d_E.mat', participant)])
    
    % Store EEG and Resp into dat and resp
    dat         = eval(sprintf('data_EEG_Suj%02d_E', participant));
    resp        = eval(sprintf('Resp_phase_Suj%02d_E', participant));

    % Clear memory
    clear(sprintf('data_EEG_Suj%02d_E', participant), sprintf('Resp_phase_Suj%02d_E', participant));
    
    % Keep the 4 s before press
    dat.trial   = cellfun(@(x) x(1:2049),dat.trial,'UniformOutput',false);
    concat_data = cell2mat(dat.trial);

    % Bin according to Resp phase
    concat_resp = reshape(resp(:,1:2049)', [1 size(resp,1)*size(resp(:,1:2049),2)]);
    
    bin_ID      = 0;
    mean_amp    = nan(1,6);
    for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
        bin_ID  = bin_ID + 1;
        I       = find(concat_resp >= phase_bin & concat_resp < phase_bin + pi/3);
        mean_amp(bin_ID) = mean(concat_data(I),'all');
    end

    % Get 1000 permutations
    perm_amp        = nan(1000,6);
    perm_omnibus    = nan(1000,1);
    
    for i = 1:1000
        temp_resp = circshift(concat_resp, randi(size(concat_resp,2)));
        bin_ID      = 0;

        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find(temp_resp >= phase_bin & temp_resp < phase_bin + pi/3);
            perm_amp(i,bin_ID) = mean(concat_data(I),'all');
        end
        [~,perm_omnibus(i)] = circ_otest(temp_resp(2049:2049:end));
    end
    [~,m] = circ_otest(resp(:,2049));
    
    % Save files
    originaldata_originalbinning(participant).mean_amp           = mean_amp;
    originaldata_originalbinning(participant).perm_amp           = perm_amp;
    originaldata_originalbinning(participant).omnibus            = m;
    originaldata_originalbinning(participant).perm_omnibus       = perm_omnibus;

    % With trial level normalization (using Park et al.'s 2020 normalization by sum division)
    temp_mean_amp       = nan(length(dat.trial),6);
    mean_amp_sumnorm    = nan(length(dat.trial),6);

    for trl = 1:length(dat.trial)
        bin_ID = 0;
        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find( resp(trl,1:2049) >= phase_bin & resp(trl,1:2049) < phase_bin + pi/3);
            temp_mean_amp(trl, bin_ID) = mean(dat.trial{trl}(I),'all');
        end
        if any(isnan(temp_mean_amp(trl, :))) % remove trials with any bins with a nan
            temp_mean_amp(trl, :)    = nan(1,6);
        else
            mean_amp_sumnorm(trl, :)    = temp_mean_amp(trl, :) / sum(temp_mean_amp(trl, :));
        end
    end

    originaldata_originalbinning(participant).mean_amp_sumnorm        = mean_amp_sumnorm; 
end

fprintf('Original data binning done. Starting permutation statistics ... \n')

%% Permutation statistics - Original Data

for participant = 1:length(originaldata_originalbinning)  
    temp_amp = originaldata_originalbinning(participant).mean_amp + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
    n = length(temp_amp);
    p = temp_amp / sum(temp_amp);
    H = -sum(p.*log(p));
    originaldata_originalbinning(participant).MI = (log(n) - H) / log(n);
    
    originaldata_originalbinning(participant).MI_perm = nan(1, 1000);

    % Permutations
    for i = 1:1000
        temp_amp = originaldata_originalbinning(participant).perm_amp(i,:) + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
        n = length(temp_amp);
        p = temp_amp / sum(temp_amp);
        H = -sum(p.*log(p));
        originaldata_originalbinning(participant).MI_perm(i) = (log(n) - H) / log(n);
    end

    % Individual significance
    originaldata_originalbinning(participant).pval_MI = sum(originaldata_originalbinning(participant).MI_perm > originaldata_originalbinning(participant).MI) / 1000; 
end

sum_MI_perm_vals    = sum(reshape([originaldata_originalbinning.MI_perm], [], length(originaldata_originalbinning))',1);
sum_MI_vals         = sum(reshape([originaldata_originalbinning.MI], [], length(originaldata_originalbinning))',1);
originaldata_originalbinning(1).pval_MI_total = (sum(sum_MI_perm_vals > sum_MI_vals) + 1) / 1001;

summed_perm_omni = sum(reshape([originaldata_originalbinning.perm_omnibus], [],length(originaldata_originalbinning)),2);
summed_omni      = sum([originaldata_originalbinning.omnibus]);
originaldata_originalbinning(1).pval_omni_total = (sum(summed_perm_omni<summed_omni) + 1) / 1001;

% Outlier rejection 3*SD method
perm_data       = reshape([originaldata_originalbinning.MI_perm], [], length(originaldata_originalbinning))';
orig_data       = reshape([originaldata_originalbinning.MI], [], length(originaldata_originalbinning))';
reject_thres    = 3*std(orig_data);
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
originaldata_originalbinning(1).pval_MI_total_outlier_3sd = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

% Outlier rejection median-IQR method
outlier_iqr     = prctile(orig_data,75)-prctile(orig_data,25);
reject_thres    = prctile(orig_data,75)+1.5*outlier_iqr;
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
originaldata_originalbinning(1).pval_MI_total_outlier_iqr = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

save('results\original_binning\originaldata_originalbinning.mat', 'originaldata_originalbinning')
fprintf('Original data binning done and saved \n')

%% ----------------------------- New Data ------------------------------ %%
clear
files_new                       = dir('data\new\processed\P*_data_cube.mat');
newdata_originalbinning         = struct();

for participant = 1:length(files_new)
    % Load data
    load([files_new(participant).folder '\' sprintf('P%02d_data_cube.mat', participant)])
    load([files_new(participant).folder '\' sprintf('P%02d_resp_segment.mat',participant)])

    % Store EEG and Resp into dat and resp.
    dat     = squeeze(data_cube(1,:,:));
    resp    = squeeze(data_cube(2,:,:));

    % Clear memory
    clear('data_cube');
      
    % Preprocess and unwrap eeg and resp for binning (taking 4s before
    % press)
    concat_data = reshape(dat(1:2049,:), [1 size(dat,2)*size(dat(1:2049,:),1)] );
    concat_resp = reshape(resp(1:2049,:), [1 size(resp,2)*size(resp(1:2049,:),1)]); 
    
    % Add or subtract pi to account for .eeg save sign inversion 
    concat_resp_tp = concat_resp;
    concat_resp(concat_resp_tp > 0) =  concat_resp_tp(concat_resp_tp > 0) - pi;
    concat_resp(concat_resp_tp < 0) =  concat_resp_tp(concat_resp_tp < 0) + pi;

    resp_tp = resp;
    for trl_ct = 1:size(resp,2)
        resp(resp_tp(:,trl_ct) > 0,trl_ct) =  resp_tp(resp_tp(:,trl_ct) > 0,trl_ct) - pi;
        resp(resp_tp(:,trl_ct) < 0,trl_ct) =  resp_tp(resp_tp(:,trl_ct) < 0,trl_ct) + pi;
    end
    

    bin_ID = 0;
    mean_amp = nan(1,6);
    for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
        bin_ID = bin_ID + 1;
        I = find(concat_resp >= phase_bin & concat_resp < phase_bin + pi/3);
        mean_amp(bin_ID) = mean(concat_data(I),'all');
    end

    % Get permutations
    perm_amp = nan(1000,6);
    
    for i = 1:1000
        temp_resp = circshift(concat_resp, randi(size(concat_resp,2)));
        temp_segment = circshift(resp_segment, randi(size(resp_segment,2)));
        bin_ID = 0;

        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find(temp_resp >= phase_bin & temp_resp < phase_bin + pi/3);
            perm_amp(i,bin_ID) = mean(concat_data(I),'all');
        end
        [~,perm_omnibus(i)] = circ_otest(temp_segment(press_positions));
    end
    [~,m] = circ_otest(resp_segment(press_positions));

    % Save files
    newdata_originalbinning(participant).mean_amp           = mean_amp;                % Used for Bayesian test and MI test
    newdata_originalbinning(participant).perm_amp           = perm_amp;                % Used for MI test
    newdata_originalbinning(participant).omnibus            = m;
    newdata_originalbinning(participant).perm_omnibus       = perm_omnibus;

    temp_mean_amp       = nan(size(dat,2),6);
    mean_amp_sumnorm    = nan(size(dat,2),6);

    for trl = 1:size(dat,2)
        bin_ID = 0;
        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find( resp(1:2049, trl) >= phase_bin & resp(1:2049, trl) < phase_bin + pi/3);
            temp_mean_amp(trl, bin_ID) = mean(dat(I, trl),'all');
        end
        if any(isnan(temp_mean_amp(trl, :))) % remove trials with any bins with a nan
            temp_mean_amp(trl, :)    = nan(1,6);
        else
            mean_amp_sumnorm(trl, :)    = temp_mean_amp(trl, :) / sum(temp_mean_amp(trl, :)); % Park et al. sum normalization
        end
    end

    newdata_originalbinning(participant).mean_amp_sumnorm        = mean_amp_sumnorm;        % Use for appendix reproduce figure

end

fprintf('New data binning done. Starting permutation statistics ... \n')

%% Permutation statistics - New Data

for participant = 1:length(newdata_originalbinning)  
    temp_amp = newdata_originalbinning(participant).mean_amp + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
    n = length(temp_amp);
    p = temp_amp / sum(temp_amp);
    H = -sum(p.*log(p));
    newdata_originalbinning(participant).MI = (log(n) - H) / log(n);
    
    newdata_originalbinning(participant).MI_perm = nan(1, 1000);

    % Permutations
    for i = 1:1000
        temp_amp = newdata_originalbinning(participant).perm_amp(i,:) + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
        n = length(temp_amp);
        p = temp_amp / sum(temp_amp);
        H = -sum(p.*log(p));
        newdata_originalbinning(participant).MI_perm(i) = (log(n) - H) / log(n);
    end

    % Individual significance
    newdata_originalbinning(participant).pval_MI = sum(newdata_originalbinning(participant).MI_perm > newdata_originalbinning(participant).MI) / 1000;
end

sum_MI_perm_vals    = sum(reshape([newdata_originalbinning.MI_perm], [], length(newdata_originalbinning))',1);
sum_MI_vals         = sum(reshape([newdata_originalbinning.MI], [], length(newdata_originalbinning))',1);
newdata_originalbinning(1).pval_MI_total = (sum(sum_MI_perm_vals > sum_MI_vals) + 1) / 1001;

summed_perm_omni = sum(reshape([newdata_originalbinning.perm_omnibus], [],length(newdata_originalbinning)),2);
summed_omni      = sum([newdata_originalbinning.omnibus]);
newdata_originalbinning(1).pval_omni_total = (sum(summed_perm_omni<summed_omni) + 1) / 1001;

% Outlier rejection 3*SD method
perm_data       = reshape([newdata_originalbinning.MI_perm], [], length(newdata_originalbinning))';
orig_data       = reshape([newdata_originalbinning.MI], [], length(newdata_originalbinning))';
reject_thres    = 3*std(orig_data);
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
newdata_originalbinning(1).pval_MI_total_outlier_3sd = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

% Outlier rejection median-IQR method
outlier_iqr     = prctile(orig_data,75)-prctile(orig_data,25);
reject_thres    = prctile(orig_data,75)+1.5*outlier_iqr;
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
newdata_originalbinning(1).pval_MI_total_outlier_iqr = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

save('results\original_binning\newdata_originalbinning.mat', 'newdata_originalbinning')
fprintf('New data binning done and saved \n')


%% -------------------------- Simulated Data --------------------------- %%
clear
files_simulated                  = dir('data\simulated\processed\SIMUL*_cube.mat');
simulateddata_originalbinning    = struct();

for participant = 1:length(files_simulated)
    % Load data
    load([files_simulated(participant).folder '\' sprintf('SIMUL%02d_cube.mat', participant)])

    % Extract dat and resp
    dat     = squeeze(data_cube(1,:,:));
    resp    = squeeze(data_cube(2,:,:));

    % Clear memory
    clear('data_cube');
      
    % Keep 4 s before press
    concat_data = reshape(dat(1:2048,:), [1 size(dat,2)*size(dat(1:2048,:),1)] );

    % Bin according to Resp phase
    concat_resp = reshape(resp(1:2048,:), [1 size(resp,2)*size(resp(1:2048,:),1)]);
    
    bin_ID = 0;
    mean_amp = nan(1,6);
    for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
        bin_ID = bin_ID + 1;
        I = find(concat_resp >= phase_bin & concat_resp < phase_bin + pi/3);
        mean_amp(bin_ID) = mean(concat_data(I),'all');
    end

    % Get permutations
    perm_amp = nan(1000,6);
    
    for i = 1:1000
        temp_resp = circshift(concat_resp, randi(size(concat_resp,2)));
        bin_ID = 0;

        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find(temp_resp >= phase_bin & temp_resp < phase_bin + pi/3);
            perm_amp(i,bin_ID) = mean(concat_data(I),'all');
        end
        [~,perm_omnibus(i)] = circ_otest(temp_resp(2048:2048:end));
    end
    [~,m] = circ_otest(resp(2048,:));
    
    % Save files
    simulateddata_originalbinning(participant).mean_amp           = mean_amp;                % Used for Bayesian test and MI test
    simulateddata_originalbinning(participant).perm_amp           = perm_amp;                % Used for MI test
    simulateddata_originalbinning(participant).omnibus            = m;
    simulateddata_originalbinning(participant).perm_omnibus       = perm_omnibus;
    
    % With trial level normalization (using Park et al.'s sum)
    temp_mean_amp       = nan(size(dat,2),6);
    mean_amp_sumnorm    = nan(size(dat,2),6);

    for trl = 1:size(dat,2)
        bin_ID = 0;
        for phase_bin = [-pi -pi*(2/3) -pi/3 0 pi/3 2*pi/3]
            bin_ID = bin_ID + 1;
            I = find( resp(1:2048, trl) >= phase_bin & resp(1:2048, trl) < phase_bin + pi/3);
            temp_mean_amp(trl, bin_ID) = mean(dat(I, trl),'all');
        end
        if any(isnan(temp_mean_amp(trl, :))) % remove trials with any bins with a nan
            temp_mean_amp(trl, :)    = nan(1,6);
        else
            mean_amp_sumnorm(trl, :)    = temp_mean_amp(trl, :) / sum(temp_mean_amp(trl, :)); % Park et al. sum normalization
        end
    end

    simulateddata_originalbinning(participant).mean_amp_sumnorm        = mean_amp_sumnorm;        % Use for appendix reproduce figure
end

fprintf('Simulated data binning done. Starting permutation statistics ... \n')

%% Permutation statistics - Simulated Data

for participant = 1:length(simulateddata_originalbinning)  
    temp_amp = simulateddata_originalbinning(participant).mean_amp + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
    n = length(temp_amp);
    p = temp_amp / sum(temp_amp);
    H = -sum(p.*log(p));
    simulateddata_originalbinning(participant).MI = (log(n) - H) / log(n);
    
    simulateddata_originalbinning(participant).MI_perm = nan(1, 1000);

    % Permutations
    for i = 1:1000
        temp_amp = simulateddata_originalbinning(participant).perm_amp(i,:) + 500; % Fix from Park et al. 2020 to avoid negative in the log(p) (see their code for justification)
        n = length(temp_amp);
        p = temp_amp / sum(temp_amp);
        H = -sum(p.*log(p));
        simulateddata_originalbinning(participant).MI_perm(i) = (log(n) - H) / log(n);
    end

    % Individual significance
    simulateddata_originalbinning(participant).pval_MI = sum(simulateddata_originalbinning(participant).MI_perm > simulateddata_originalbinning(participant).MI) / 1000;
end

sum_MI_perm_vals    = sum(reshape([simulateddata_originalbinning.MI_perm], [], length(simulateddata_originalbinning))',1);
sum_MI_vals         = sum(reshape([simulateddata_originalbinning.MI], [], length(simulateddata_originalbinning))',1);
simulateddata_originalbinning(1).pval_MI_total = (sum(sum_MI_perm_vals > sum_MI_vals) + 1) / 1001;

% Outlier rejection
perm_data       = reshape([simulateddata_originalbinning.MI_perm], [], length(simulateddata_originalbinning))';
orig_data       = reshape([simulateddata_originalbinning.MI], [], length(simulateddata_originalbinning))';
reject_thres    = 3*std(orig_data);
rej_idx         = orig_data>reject_thres;
sum_MI_perm_vals_out    = sum(perm_data(~rej_idx,:),1);
sum_MI_vals_out         = sum(orig_data(~rej_idx),1);
simulateddata_originalbinning(1).pval_MI_total_outlier = (sum(sum_MI_perm_vals_out > sum_MI_vals_out) + 1) / 1001;

save('results\original_binning\simulateddata_originalbinning.mat', 'simulateddata_originalbinning')
fprintf('Simulated data binning done and saved \n')