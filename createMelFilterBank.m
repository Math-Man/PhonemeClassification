function [filterBank] = createMelFilterBank(frames, sampleRate, lowerBound, upperBound, filterCount, nfft)


lowFreqMel = (2595 * log10(1 + (lowerBound / 2) / 700));
highFreqMel = (2595 * log10(1 + (upperBound / 2) / 700));

melPoints = linspace(lowFreqMel, highFreqMel, filterCount+2);

hertzPoints = (700 * (10.^(melPoints/2595) - 1));

%Round down the values to fit within the resolution band width
melScaledBins = floor( (nfft + 1) * hertzPoints / sampleRate);

fBank = zeros(filterCount, floor(nfft/2 + 1));

hold on
for m = (2:(filterCount+1))
    
    fprev = melScaledBins(m-1);
    fcurrent = melScaledBins(m);
    fnext = melScaledBins(m+1);
    
    for k = fprev:fcurrent
        fBank(m - 1, k) = (k - melScaledBins(m-1)) / (melScaledBins(m) -  melScaledBins(m - 1));
    end
    
    for k = fcurrent:fnext
        fBank(m - 1, k) = (melScaledBins(m+1) - k) / (melScaledBins(m+1) -  melScaledBins(m));
    end
    plot(fBank(m-1,:));
end
hold off


filterBank = fBank;

end
% 
% 
% 
% 
% 
% %Periodgram
% frameFFT = [];
% for frame = 1:length(frames) 
%     frameFFT(frame,:)= (fft(frames(frame),nfft)).^2; 
% end
% 
% %mel scale
% melScale = [];
% count = 1;
% for i = lowerBound:( floor((upperBound-lowerBound)/(filterCount+1)) ):upperBound%Creates first M+2 filters
%     melScale(count) = 1125 * log(1 + i/700);
%     count = count + 1;
% end
% 
% 
% %inverse mel scale
% melScaleInverse = [];
% for i = 1:length(melScale)
%     melScaleInverse(i) = 700 * (exp(melScale(i)/1125) - 1);
% end
% 
% 
% 
% 
% 



%Create the filterbanks using the following formulation:
%   for filterbank Hm(k)
%   0 if k < melScaledFrequencies(m - 1)
%   (k - melScaledFrequencies(m - 1)) / (melScaledFrequencies(m) - melScaledFrequencies(m - 1)) for melScaledFrequencies(m-1) <= k <= melScaledFrequencies(m)
%   (k + melScaledFrequencies(m + 1)) / (melScaledFrequencies(m+1) - melScaledFrequencies(m)) for melScaledFrequencies(m) <= k <= melScaledFrequencies(m+1)
%   0 if k > melScaledFrequencies(m + 1)
%   This algorithm uses triang(L) for calculation
%One filter bank's peak starts the other filter



% melFilterBank = [];
% filterIndex = 1;
% filterStart = 0;
% 
% for i = 1:length(melScaledFrequencies)
%     filterShape = triang(melScaledFrequencies(i));  %Triang with length of the melscaledfreqs.
%     for k = 1:melScaledFrequencies(i)
%         
%         melFilterBank(i, filterIndex - filterStart) = filterShape(k);
%         
%         filterIndex = filterIndex + 1;
%     end
%     filterStart = ceil(melScaledFrequencies(i)/2) + filterStart;
% end
% 
% hold on
% for i = 1:length(melScaledFrequencies)
%     plot(melFilterBank(i,:));
% end
% hold off



% filterStart = 1;
% for filterIndex = 2: length(melScaledFrequencies)-1    
%     for k = filterStart : (filterStart+(melScaledFrequencies(filterIndex)))
%         
%         %Previous filter
%             previousFilter = melScaledFrequencies(filterIndex - 1)
%         
%         %CurrentFilter
%         currentFilter = melScaledFrequencies(filterIndex);
%         
%         %next filter
%             nextFilter =  melScaledFrequencies(filterIndex + 1);
% 
%         
%         peakFound = 0;
%         
%         if (k < previousFilter) 
%             melFilters(k) = 0;
%             
%         elseif (previousFilter <= k && k <= currentFilter )
%             
%             melFilters(k) = (k - previousFilter)/(currentFilter - previousFilter);
%             
%         elseif (currentFilter <= k && k <= nextFilter )
%             
%             if(peakFound == 0)
%                 peakFound = 1;
%                 filterStart = k;
%             end
%             
%             melFilters(k) = (k + nextFilter)/(nextFilter - currentFilter);
%             
%         else
%             melFilters(k) = 0;
%         end
%     end
%     
%     plot(melFilters);
%     
% end

