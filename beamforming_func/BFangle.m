function out = BFangle(in, maxAprSz, Fs, c, pitch, apo, Fnum, theta)
%% BFANGLE Beamforms with angular steering
% function out = BFangle(in, maxAprSz, Fs, c, pitch, apo, Fnum, theta)
% The following script takes the pre-beamform data and creates a beamformed
% image by applying parallel beamformers to the same data. 
%
% input:
%   in:             the first dimension should be scan lines and the second
%                   dimension should be the channel samples
%   maxAprSz:       maximum size of the aperture for beamforming
%   Fs              Acquisition sampling frequency 80*1e6, 40MHz, ... 10*1e6
%   c               speed of sound [m/s]
%   pitch           distance between elements [m]
% output:
%   out: beamformed rf data
%
% Example:
%           [hdr, RF] = readDAQ(path, chanls, frameN , reRoute);
%           bRF = BF(RF, 64, 40*1e6, 2, 1540, .3048E-3, hanning(64));
%s
% see also
%       readDAQ.m
%
% Author: Reza Zahiri Azar
% Copyright Ultrasonix Medical Corporation - Dec 2011


% The first dimension has to be lines and second dimension has to be samples
if (size(in,1) > size(in,2))
    in = in';
    transposed = true;
else
    transposed = false;
end

nl = size(in,1);   % number of lines
ns = size(in,2);   % number of samples

sampleSpacing = c/Fs/2;     % spacing between samples in the axial direction in milimeter for 40MHz

offset = 0.00;                 % distance between first sample and transducer element [m]
nSampleOffset = round(offset / sampleSpacing);    % offset in samples: 1st sample is actually 1+nSampleOffset sample

out = zeros(nl,ns);

%hlfAprSz = floor(maxAprSz / 2);  
%x = -hlfAprSz: hlfAprSz;    % aperture indices
%fulAprSz = 2*hlfAprSz + 1;


parfor i = 1:nl    % for each line/channel
    % disp(['scanline #', num2str(i)]);
    
    for j=1:ns  % find the value for each samples

        % find the sample location
        Z =  (j + nSampleOffset) * sampleSpacing ; 

        % calculate the aperture based on the F number
        a = Z/(2*Fnum);
        hlfAprSz = floor(a / pitch);
        % chk to make sure we do not go beyound max apr
        if (hlfAprSz > maxAprSz/2)
            hlfAprSz = floor(maxAprSz / 2);
        end
        x = -hlfAprSz: hlfAprSz;    % aperture indices
        fulAprSz = 2*hlfAprSz + 1;

        if strcmp(apo, 'rect')
            win = ones(1, fulAprSz);
        elseif strcmp(apo,'bh')
            win = blackmanharris(fulAprSz); win = win(:)';
        elseif strcmp(apo, 'nullRECT')
            win = [ones(1,floor(fulAprSz/2)) 0 -1*ones(1,floor(fulAprSz/2))];
        elseif strcmp(apo, 'nullRECTl')
            win = [ones(1,floor(fulAprSz/2)) 0 -1*ones(1,floor(fulAprSz/2))] + apod_dc;
        elseif strcmp(apo, 'nullRECTr')
            win = fliplr([ones(1,floor(fulAprSz/2)) 0 -1*ones(1,floor(fulAprSz/2))] + apod_dc);
        elseif strcmp(apo, 'nullRECTbh')
            win = blackmanharris(fulAprSz)'.*[ones(1,floor(fulAprSz/2)) 0 -1*ones(1,floor(fulAprSz/2))];
        elseif strcmp(apo, 'nullRECTlbh')
            win = blackmanharris(fulAprSz)'.*[ones(1,floor(fulAprSz/2)) 0 -1*ones(1,floor(fulAprSz/2))] + apod_dc;
        elseif strcmp(apo, 'nullRECTrbh')
            win = blackmanharris(fulAprSz)'.*fliplr([ones(1,floor(fulAprSz/2)) 0 -1*ones(1,floor(fulAprSz/2))] + apod_dc);s
        end

%         if (i==83&&j==909)
%             keyboard;
%         end
        
        delays = calc_delays(c,Fs,pitch,i,j,theta,fulAprSz);
        
        chnls = zeros(1,fulAprSz);  % will be filled with proper signal from each channel
        
        cntr = i;       % center of aperture
        apr = cntr + x; % absolute index locations

        % find the corresponding value from each channel
        for k = 1:fulAprSz  
            
            % clip to data boundaries
            %   lateral:
            if apr(k)<1, continue, end;
            if apr(k)>nl, continue, end;
            
            %   axial:
            if delays(k)<1, continue, end;
            if delays(k)>ns, continue, end;
            
            chnls(k) = in(apr(k), delays(k)); 
            
        end;
        % apodization : ideally has to be a function of depth
        chnls = win .* chnls;
        
        % beamforming
        % if wanted digital TGC can be applied here e.g TGC = 1 + 2 * (j/ns);
        out(i,j) = sum( chnls );
        
    end;
end

if (transposed),    out = out';
end
