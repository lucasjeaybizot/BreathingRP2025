function [resp_EMG] = GetEMG()
%GETEMG Extracts respiration phase at EMG onset for the new dataset.
%
% This function processes EEG data files for 17 participants. For each participant,
% it:
%  - Loads raw EEG data,
%  - Downsamples the signal to 512 Hz,
%  - Filters and extracts the respiration phase using the Hilbert transform,
%  - Filters and rectifies EMG data from the right hand,
%  - Identifies EMG onset using a percentile-based method near known stimulus events,
%  - Extracts the respiration phase at the adjusted EMG onset,
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
    
        %% Bandpass filter & rectify EMG data
        cfg             = [];
        cfg.channel     = find(strcmp(data_rs.label, 'R Hand'));
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [20 200];
        data_emg_bp     = ft_preprocessing(cfg, data_rs);
        data_emg_bp.trial{1,1} = lowpass(abs(data_emg_bp.trial{1,1}),10,512); 
    
        %% Identify EMG onset
        ev              = ft_read_event(filename);
        ev              = ev(strcmp({ev.type}, 'Stimulus'));
        trial_idx       = [find(strcmp({ev.value}, 'S 20'),50,'first') ...
                           find(strcmp({ev.value}, 'S 30'),50,'first')];
    
        press_samples   = round([ev(trial_idx).sample]/(2500/512));
    
        for trl = 1:length(press_samples)
            first500ms = data_emg_bp.trial{1,1}((press_samples(trl)-255):(press_samples(trl)));
            jittered_emg = 256-find(first500ms>prctile(data_emg_bp.trial{1,1}(max([(press_samples(trl)-1023) 1]):(press_samples(trl))),97.5),1,'first');
            if isempty(jittered_emg)
                jittered_emg=0;
            end
            resp_EMG(participant).resp_phases(trl) = data_resp_bp.trial{1,1}(press_samples(trl)-jittered_emg);
        end
    
    end

    save('results\figure_S3\resp_EMG.mat', 'resp_EMG')

end

