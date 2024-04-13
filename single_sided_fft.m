function [fft_positive, frequencies] = single_sided_fft(signal, sample_rate)
    % calculate single-sided amplitude spectrum
    % Inputs:
    %   signal: The input signal
    %   sample_rate: The sample rate of the signal
    % Outputs:
    %   fft_positive: The single-sided FFT amplitude spectrum
    %   frequencies: The corresponding frequencies
    
    % get signal length
    signal_length = numel(signal);
  
    % apply fft
    fft_signal = fft(signal);
    % scale signal by signal length 
    fft_scaled = fft_signal/signal_length;
    % take magnitude
    fft_magnitude = abs(fft_scaled);
    % cut off negative frequencies
    fft_positive = fft_magnitude(1:signal_length/2+1);
    % multiply positive frequencies by 2
    fft_positive(2:end-1) = 2*fft_positive(2:end-1);
    
    % define frequencies
    frequencies = sample_rate/signal_length * (0:(signal_length/2));

end

