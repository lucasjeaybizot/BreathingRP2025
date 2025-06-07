clear

data_cube = nan(2,2048,100); % initiate empty 3D matrix

for j = 1:100
    % Set the pseudo-EEG signal
    data_cube(1,:,j) = [zeros(1,1024) linspace(0,-4,1024)]; % set the exact same pseudo-RP for all trials

    % Set the phase signal per trial
    % At each if-statement the last item added inside the circshift sets the phase locking of the selected trials
    % 2048 sets the sampling rate
    
    if j<50  % select 50 percent of trials
            data_cube(2,:,j) = circshift(linspace(-pi,pi,2048), round(2048/12) + 0);
    elseif j>=50 & j<60
            data_cube(2,:,j) = circshift(linspace(-pi,pi,2048), round(2048/12) + round(2048/6));
    elseif j>=60 & j<70
            data_cube(2,:,j) = circshift(linspace(-pi,pi,2048), round(2048/12) + round(2*2048/6));
    elseif j>=70 & j<80
            data_cube(2,:,j) = circshift(linspace(-pi,pi,2048), round(2048/12) + round(3*2048/6));
    elseif j>=80 & j<90
            data_cube(2,:,j) = circshift(linspace(-pi,pi,2048), round(2048/12) + round(4*2048/6));
    elseif j>=90
            data_cube(2,:,j) = circshift(linspace(-pi,pi,2048), round(2048/12) + round(5*2048/6));
    end
end

% Copy this simulated data 50 times for 50 pseudo participants
for i = 1:50
    save(sprintf('data\\simulated\\processed\\SIMUL%02d_cube.mat', i), 'data_cube');
end