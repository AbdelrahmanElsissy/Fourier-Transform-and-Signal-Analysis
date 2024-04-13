function ramp_filter = create_ramp_filter(ramp_length, signal_length)
    % creates symmetric ramp-filter based on hanning window filter
    % Inputs:
    %   ramp_length: Length of the ramping for start and ending
    %   signal_length: Length of the input signal
    % Outputs:
    %   ramp_filter: Ramp-filter which is a matrix of ones except the start
    %   and end with the defined ramping length

    % create hann filter
    hann_filter = hann(2*ramp_length);
    % split hann filter into two sections
    ramp = hann_filter(1:ramp_length);
    deramp = hann_filter(ramp_length+1:end);

    % create ramp filter with signal_length
    ramp_filter = ones([signal_length, 1]);
    % replace beginning and end with the ramping
    ramp_filter(1:ramp_length) = ramp;
    ramp_filter(signal_length-ramp_length+1:end) = deramp;

end