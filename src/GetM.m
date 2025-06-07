function [resp_M] = GetM()
%GETM Extracts respiration phase for M-trials only for the new dataset.
%
% This function processes EEG data files for 17 participants. For each participant,
% it:
%  - Loads raw EEG data,
%  - Downsamples the signal to 512 Hz,
%  - Filters and extracts the respiration phase using the Hilbert transform,
%  - Extracts the respiration phase at button press only for M trials,
%  - Stores results in a struct and saves them to results\figure_S3 folder.

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
    
        %% Bandpass filter respiration (0.2-0.8 Hz)
        cfg             = [];
        cfg.channel     = find(strcmp(data_rs.label, 'Resp'));
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [0.2 0.8];
        cfg.bpfiltord   = 3;
        data_resp_bp    = ft_preprocessing(cfg, data_rs);
    
        %% Extract respiration phase using Hilbert transform
        data_resp_bp.trial{1,1} = angle(hilbert(data_resp_bp.trial{1,1})); 
    
        % Account for sign inversion when saving .eeg data
        dat_tp          = data_resp_bp.trial{1,1};
        data_resp_bp.trial{1,1}(dat_tp>0) = dat_tp(dat_tp>0) - pi;
        data_resp_bp.trial{1,1}(dat_tp<0) = dat_tp(dat_tp<0) + pi; 
    
        %% Identify EMG onset
        ev              = ft_read_event(filename);
        ev              = ev(strcmp({ev.type}, 'Stimulus'));
        trial_idx       = [find(strcmp({ev.value}, 'S 20'),50,'first')];
    
        press_samples   = round([ev(trial_idx).sample]/(2500/512));
    
        resp_M(participant).resp_phases = data_resp_bp.trial{1,1}(press_samples);
    
    end

    save('results\figure_S3\resp_M.mat', 'resp_M')

end

