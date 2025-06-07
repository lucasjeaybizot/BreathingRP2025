%% Preprocessing EEG and Respiration Data for new data
% This script processes EEG and respiration data for the 17 participants
% from the new dataset. The original dataset was already preprocessed. It
% performs preprocessing steps as described in Park et al. (2020). The
% Fieldtrip toolbox is required.
% Final output: a data cube of EEG and Respiration epochs and a respiration
% segment for permutation analysis.

close all
clear

%% Loop over participants
for participant = 1:17
    
    fprintf('Processing participant %02d ... \n', participant)

    %% Load new participant
    filename        = sprintf('data/new/raw/P%02d_RespRP.eeg', participant);
    cfg             = [];
    cfg.dataset     = filename;
    data            = ft_preprocessing(cfg);

    %% Downsample data to 512 Hz
    cfg             = [];
    cfg.resamplefs  = 512;
    data_rs         = ft_resampledata(cfg, data);

    %% Bandpass filter EEG (0.1-40 Hz)
    cfg             = [];
    cfg.channel     = 1:63;
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [0.1 40];
    data_eeg_bp     = ft_preprocessing(cfg, data_rs);

    %% Bandpass filter respiration (0.2-0.8 Hz)
    cfg             = [];
    cfg.channel     = find(strcmp(data_rs.label, 'Resp'));
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [0.2 0.8];
    cfg.bpfiltord   = 3;
    data_resp_bp    = ft_preprocessing(cfg, data_rs);

    %% Extract respiration phase using Hilbert transform
    data_resp_bp.trial{1,1} = angle(hilbert(data_resp_bp.trial{1,1})); 

    % Account for sign inversion in .eeg saving
    dat_tp = data_resp_bp.trial{1,1};
    data_resp_bp.trial{1,1}(dat_tp>0) = dat_tp(dat_tp>0) - pi;
    data_resp_bp.trial{1,1}(dat_tp<0) = dat_tp(dat_tp<0) + pi;

    %% Rereference EEG to common average
    cfg             = [];
    cfg.reref       = 'yes';
    cfg.refchannel  = 'all';
    cfg.implicitref = 'Oz';
    data_eeg_reref  = ft_preprocessing(cfg, data_eeg_bp);

    %% Epoch around button press [-4s, +1s]
    ev              = ft_read_event(filename);
    ev              = ev(strcmp({ev.type}, 'Stimulus'));
    cfg             = [];
    trial_idx       = [find(strcmp({ev.value}, 'S 20'),50,'first') ...
                       find(strcmp({ev.value}, 'S 30'),50,'first')];
    
    % Construct trial definition matrix (trl)
    fs              = data_eeg_reref.fsample;
    press_samples   = round([ev(trial_idx).sample]/(2500/512));
    cfg.trl         = [-4*fs + press_samples; ...
                       +1*fs + press_samples; ...
                       -4*fs*ones(1,length(trial_idx))]';

    % Ensure all trials are within bounds
    if cfg.trl(1,1) < 0
        cfg.trl(1,:) = [];
        trial_idx(1) = [];
        press_samples(1) = [];
    end

    % Epoch EEG and respiration
    data_eeg_epoch  = ft_redefinetrial(cfg, data_eeg_reref);
    data_resp_epoch = ft_redefinetrial(cfg, data_resp_bp);

    % Store respiration and press position for permutation tests
    segment_resp_idx= [min(cfg.trl(:,1)) max(cfg.trl(:,2))];
    press_positions = press_samples - min(cfg.trl(:,1));
    resp_segment    = data_resp_bp.trial{1,1}(1,segment_resp_idx(1):segment_resp_idx(2));
    
    %% Reject trials with abnormal amplitude (>3 SD)
    chan_of_interest= ismember(data_eeg_epoch.label, {'Cz', 'FCz','Fz','AFz'});
    to_reject       = zeros(size(trial_idx));
    trl_mean_amp    = zeros(4, size(trial_idx,2));

    for i = 1:size(trial_idx,2)
        trl_mean_amp(:, i)  = mean(data_eeg_epoch.trial{1,i}(chan_of_interest,:), 2);
    end

    for i = 1:size(trial_idx,2)
        to_reject(i) = any(abs(trl_mean_amp(:,i)) > 3*std(trl_mean_amp,[],2), 'all');
    end
    
    press_positions = press_positions(~to_reject);

    %% Extract RP trials from the most negative channel among Cz, FCz, Fz and AFz
    chan_idx        = find(chan_of_interest);
    RP_Fz           = nan(1, 2561, sum(to_reject==0));
    RP_Cz           = nan(1, 2561, sum(to_reject==0));
    RP_AFz          = nan(1, 2561, sum(to_reject==0));
    RP_FCz          = nan(1, 2561, sum(to_reject==0));

    trl_I = 0;
    for i = 1:length(to_reject)
        if to_reject(i) == 0
            trl_I = trl_I + 1;
            RP_Fz(1,:,trl_I)    = data_eeg_epoch.trial{1,i}(chan_idx(1),:);
            RP_Cz(1,:,trl_I)    = data_eeg_epoch.trial{1,i}(chan_idx(2),:);
            RP_AFz(1,:,trl_I)   = data_eeg_epoch.trial{1,i}(chan_idx(3),:);
            RP_FCz(1,:,trl_I)   = data_eeg_epoch.trial{1,i}(chan_idx(4),:);
        end
    end

    % Choose channel with most negative RP amplitude (-2s to 0s = 1024:2048)
    [~,I] = min([nanmean(RP_Fz(1,1024:2048,:),[2 3]), ...
                 nanmean(RP_Cz(1,1024:2048,:),[2 3]), ...
                 nanmean(RP_AFz(1,1024:2048,:),[2 3]), ...
                 nanmean(RP_FCz(1,1024:2048,:),[2 3])]);
    
    %% Build data cube (EEG + Respiration) for retained trials
    data_cube       = nan(2,2561,sum(to_reject==0));
    trl_I = 0;
    for i = 1:length(to_reject)
        if to_reject(i) == 0
            trl_I = trl_I + 1;
            data_cube(1,:,trl_I) = data_eeg_epoch.trial{1,i}(chan_idx(I),:);
            data_cube(2,:,trl_I) = data_resp_epoch.trial{1,i};
        end
    end

    %% Save results
    save(sprintf('data/new/processed/P%02d_data_cube.mat', participant), 'data_cube');
    save(sprintf('data/new/processed/P%02d_resp_segment.mat',participant),'resp_segment','press_positions');

end
