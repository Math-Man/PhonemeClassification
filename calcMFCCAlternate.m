function [FoundCoeffs, FoundCoeffsx, frameCount, framePhonemes, countedSamples, PhonemeIndex] = calcMFCCAlternate(recording, fs, chunkSizeSeconds, frameShiftSeconds, filterbank, nfft, phonemeStarts, phonemeEnds, phonemeMatches, countedSamples, PhonemeIndex, phoneLookup)

frameLength = chunkSizeSeconds*fs;
recordingLength = length(recording);
frameShiftLength = frameShiftSeconds*fs;
frameShiftCount = ceil(frameShiftLength);

frameCount = floor(recordingLength/(frameShiftCount))-1;


%%
%Phoneme Extraction
frames = [];
framePhonemes = [];
for  frame=1:frameCount

    frameStart = (frame-1)*((frameShiftCount)) + 1;
    frameEnd = frameStart + frameLength-1;
    
    if(frameEnd > recordingLength)
        recording = padarray(recording,[(-(recordingLength- frameEnd)) 0],0,'post');
        %frames(frame,:) = (recording(frameStart:frameEnd).*hamming(frameLength));
        %break;
    end
    
    %%
    %preEmphasis
    framef = filter(1,[1, -0.97],recording(frameStart:frameEnd));
    %%
    
    frames(frame,:) = (framef.*hamming(frameLength));
    
    
    
    %Set phonemes here
    Fstart = frameStart + countedSamples;
    Fend = frameEnd + countedSamples;
    
    if(PhonemeIndex > length(phonemeStarts))
        
        framePhonemes = [framePhonemes nan];    %No more phonemes left to check, everything else is nan
    else
        if( (Fstart >= (phonemeStarts(PhonemeIndex))*fs) && (Fend <= (phonemeEnds(PhonemeIndex))*fs) ) %Current frame bounds are within the phoneme bounds
            fou = find(contains(phoneLookup{1},phonemeMatches(PhonemeIndex)));
            framePhonemes = [framePhonemes fou(1)];
        elseif (Fend > (phonemeEnds(PhonemeIndex)*fs)) %Current ending sample of the frame is beyond the current phoneme bounds
            PhonemeIndex = PhonemeIndex+1;
            framePhonemes = [framePhonemes nan];
        else
            framePhonemes = [framePhonemes nan];
        end
        %countedSamples = countedSamples + frameLength;%(frameEnd - frameStart); %add length of this frame 
    end
    
end
countedSamples = countedSamples + recordingLength;
f = sprintf('Built frames: %d / %d | Empty Samples Added: %d', frame+1, frameCount, (recordingLength- frameEnd));
disp(f);
%%

% %MFCC FEATURE EXTRACTION
% %  http://practicalcryptography.com/miscellaneous/machine-learning/%guide-mel-frequency-cepstral-coefficients-mfccs/
% 
%     dctm = @( N, M )( sqrt(2.0/M) * cos( repmat([0:N-1].',1,M).* repmat(pi*([1:M]-0.5)/M,N,1) ) );
%     % Cepstral lifter routine (see Eq. (5.12) on p.75 of [1])
%     ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );
% 
%         %% FEATURE EXTRACTION 
%     K = nfft/2;
%     L = 22;            % cepstral sine lifter parameter
%     MAG = abs( fft(frames,nfft,1) ); 
%     
%     % Triangular filterbank with uniformly spaced filters on mel scale
%     %H = trifbank( M, K, R, fs, hz2mel, mel2hz ); % size of H is M x K 
%     H = filterbank;
%     
%     
%     POW = (log10(((1/nfft) * (sum(MAG).^2))));
%   
%     
%     % Filterbank application to unique part of the magnitude spectrum
%     FBE = H * MAG(1:K,:); % FBE( FBE<1.0 ) = 1.0; % apply mel floor
%     % DCT matrix computation
%     DCT = dctm(26, 26);
%     % Conversion of logFBEs to cepstral coefficients through DCT
%     CC =  DCT * log( FBE );
%     % Cepstral lifter computation
%     lifter = ceplifter( 26, L );
%     % Cepstral liftering gives liftered cepstral coefficients
%     CC = diag( lifter ) * CC; % ~ HTK's MFCCs
%     FoundCoeffs = [POW; CC(2:end,:)];
%%
%     



%Magnitude spectrum
mag_frames = abs(fft(frames', nfft));
pow_frames = ((1/nfft) * ((mag_frames).^2));
pow_frames = pow_frames(1:length(pow_frames)/2,:);  %Crop the ffts to half coefficents






%Apply filterBank to magniutde spectrum
%filterpowers = filterbank * mag_frames(1:length(mag_frames)/2,:);%pow_frames;
filterpowers = filterbank * pow_frames;

%DCT matrix Computation
fmDCT = dctmtx(40);

%filterbankEnergy = log(filterpowers);
filterbankEnergy = log(filterpowers);




%conversion to cepstral coeff 
CC = fmDCT * filterbankEnergy;



%FoundCoeffs = dct(filterbankEnergy);%/30;

FoundCoeffs = [(5*log(sum(frames'.^2))); CC(2:end,:)]; %(1:13,:)
%FoundCoeffs = FoundCoeffs/3.5;

% meanF = mean(FoundCoeffs);
% stdF = std(FoundCoeffs);
% 
% FoundCoeffs = (FoundCoeffs - meanF)./stdF;

FoundCoeffsx = mfcc(recording, fs, 'NumCoeffs', 40, 'FFTLength', nfft, 'OverlapLength', (400-160), 'WindowLength', 400)';

figure(8376)
subplot(211)
plot(FoundCoeffs(1:13,:)')
subplot(212)
plot(FoundCoeffsx(2:14,:)')

end

