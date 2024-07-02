function delays = calc_delays(c, Fs, dx, elem_num, sample, theta, apert_size)
    
    %% Handle arguments
%     if ~exists('offset')
        offset = 0;
%     end
    
    %define aperture size
    % Slow, uncomment to restore functionality
%     if ~exist('apert_size','var')
%         apert_size = 128;
%     end
    hlfAprSz = floor(apert_size / 2);  
    aprIndices = -hlfAprSz: hlfAprSz;    % aperture indices
%     fulAprSz = 2*hlfAprSz + 1;

    %get offset in samples
    sampleSpacing = c/Fs/2;
    nSampleOffset = round(offset / sampleSpacing);
    
    %% Convention:
    
    %find sample location (m)
    Z = (sample + nSampleOffset) * sampleSpacing;
    X = (elem_num)*dx;

    %calculate delays
    X1 = (elem_num + aprIndices) * dx;
    apert_length = 127;
    if (theta >= 0)
        delays = round((Z*cosd(theta)+(X*sind(theta)) + sqrt( Z^2 + (X-X1).^2 ))/c*Fs);
    else
        delays = round((Z*cosd(theta)+(apert_length*dx-X)*sind(abs(theta)) ...
            + sqrt( Z^2 + (X-X1).^2 ))/c*Fs);
    end
    %delays = round(( Z + sqrt( Z^2 + (X-X1).^2 ) ) /c * Fs);z

end
