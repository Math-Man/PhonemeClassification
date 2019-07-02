
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2
%DATA MANAGEMENT AND CONVERSION


%Create look up table for
fileID = fopen('phonesBasic.txt','r'); %phonesBasic
ec = textscan(fileID, '%s%f');
fclose(fileID);
phoneLookup = ec;
%%

%TODO: COMBINE PHONESDETAILEDLISTS 3
%DETAILED OFFSET
%DETAILED OFFSET5 + DETAILED OFFSET
%DETAILED OFFSET6 + DETAILED OFFSET5

%Match data with labels, create X and Y arrays for model
fileID = fopen('sortedLabelBasicWithSilence02','r'); % sortedLabelsDetailed02
timesPhonemes02 = textscan(fileID, '%f%f%s');
fclose(fileID);

fileID = fopen('sortedLabelBasicWithSilence05','r'); % sortedLabelsDetailed05
timesPhonemes05 = textscan(fileID, '%f%f%s');
fclose(fileID);

fileID = fopen('sortedLabelBasicWithSilence06','r'); % sortedLabelsDetailed06
timesPhonemes06 = textscan(fileID, '%f%f%s');
fclose(fileID);

%Addup offsets
offset02_1 = timesPhonemes02{1}(end);
offset02_2 = timesPhonemes02{2}(end);

timesPhonemes05{1} = timesPhonemes05{1} + offset02_1;
timesPhonemes05{2} = timesPhonemes05{2} + offset02_2;

offset05_1 = timesPhonemes05{1}(end);
offset05_2 = timesPhonemes05{2}(end);

timesPhonemes06{1} = timesPhonemes06{1} + offset05_1;
timesPhonemes06{2} = timesPhonemes06{2} + offset05_2;

PhonemeStartTimes = [timesPhonemes02{1}(1:end); timesPhonemes05{1}(1:end); timesPhonemes06{1}(1:end)];
PhonemeEndTimes = [timesPhonemes02{2}(1:end); timesPhonemes05{2}(1:end); timesPhonemes06{2}(1:end)];
PhonemesMatching = [timesPhonemes02{3}(1:end); timesPhonemes05{3}(1:end); timesPhonemes06{3}(1:end)];


%%
%Convert phonemesLabels to their numeric counterParts  4
phonemesLabelsNumeric = zeros(1,length(PhonemesMatching));
for i = 1:length(phonemesLabelsNumeric)
     ff = find(contains(phoneLookup{1},PhonemesMatching(i)));
     phonemesLabelsNumeric(i) = ff(1);
end
'done';

%%

% 1
[recording1, fs] = audioread('FM1028_0202_063000.wav');
[recording2, fs] = audioread('FM1028_0205_063000.wav');
[recording3, fs] = audioread('FM1028_0206_063000.wav');
recording = [recording1; recording2; recording3];
clear recording1 recording2 recording3

chunkSizeSeconds = 0.025;   %Window Length in seconds 
frameShiftSeconds = 0.01;   
nfft = 8192;

filterbank = createMelFilterBankBased(fs, 10, fs, 40, nfft);


realLength = length(recording);
mfccCoeff = [];
frameMatchingPhonemes = [];
fcounter = 0;
countedSamples = 0;
PhonemeIndex = 1;
%154
nk = round(realLength/154);
%Process nk samples is one cycle, stitch
for i = 1:nk:(realLength-nk)
     
    recordingPart = recording(i:(i+nk-1));
    recordingPart = filter( [1 -0.97], 1, recordingPart );  %PreEmphasis
    [partCof, partCofx, fcount, framePhonemes, countedSamples, PhonemeIndex] = calcMFCCAlternate(recordingPart, fs, chunkSizeSeconds, frameShiftSeconds, filterbank, nfft, PhonemeStartTimes, PhonemeEndTimes, PhonemesMatching, countedSamples, PhonemeIndex, phoneLookup);
    mfccCoeff = [mfccCoeff partCof];
    frameMatchingPhonemes = [frameMatchingPhonemes framePhonemes];
    fcounter = fcounter + fcount;
    
    
    clc;
    f = sprintf('Calculating mfcc:\n%d / %d samples processed.\n%d frames built.', i-1, realLength, fcounter);
    disp(f);
    
end

%%
%Normalize Coefficents
meanF = zeros(length(mfccCoeff(:,1)),length(mfccCoeff(1,:)));
stdF = zeros(length(mfccCoeff(:,1)),length(mfccCoeff(1,:)));
normedCoeffs = zeros(length(mfccCoeff(:,1)),length(mfccCoeff(1,:)));
for i = 1:length(mfccCoeff(:,1))
    meanF(:,i) = mean(mfccCoeff(:,i));
    stdF(:,i) = std(mfccCoeff(:,i));
    
    normedCoeffs(i,:) = (mfccCoeff(i,:) -  meanF(i,:));%/ stdF(i,:);
    
end




figure(4)
hold on
for i = 1:13
    plot(normedCoeffs(i,:))
end
hold off

%%

%L1MFCC = length(mfccCoeff);

% [recordingPart, fs] = audioread('FM1028_0205_063000.wav', [i  (i + 40000)]);
% realLength = length(recording);
% %Add the 0205
% for i = 1:40000:(realLength-40000)
%     [recordingPart, fs] = audioread('FM1028_0205_063000.wav', [i  (i + 40000)]);
%     [partCof fcount] = calcMFCC(recordingPart, fs, chunkSizeSeconds, frameShiftSeconds, nfft);
%     mfccCoeff = [mfccCoeff partCof];
%     
%     fcounter = fcounter + fcount;
%     
%     clc;
%     f = sprintf('File 0205\nCalculating mfcc:\n%d / %d samples processed.\n%d frames built.', i-1+40000, realLength, fcounter);
%     disp(f);
%     
% end
% L2MFCC = length(mfccCoeff) - L1MFCC;
% 
% 
% [recordingPart, fs] = audioread('FM1028_0206_063000.wav', [i  (i + 40000)]);
% realLength = length(recording);
% %Add the 0206
% for i = 1:40000:(realLength-40000)
%     [recordingPart, fs] = audioread('FM1028_0206_063000.wav', [i  (i + 40000)]);
%     [partCof fcount] = calcMFCC(recordingPart, fs, chunkSizeSeconds, frameShiftSeconds, nfft);
%     mfccCoeff = [mfccCoeff partCof];
%     
%     fcounter = fcounter + fcount;
%     
%     clc;
%     f = sprintf('File 0206\nCalculating mfcc:\n%d / %d samples processed.\n%d frames built.', i-1+40000, realLength, fcounter);
%     disp(f);
%     
% end
% L3MFCC = length(mfccCoeff) - L1MFCC - L2MFCC;

%cof = calcMFCC(recording, fs, chunkSizeSeconds, frameShiftSeconds, nfft);


% [coeffs delta deltadelta loc] = mfcc(recording, fs, 'WindowLength', frameLength, 'OverlapLength', (frameShiftCount), 'NumCoeffs', 13, 'FFTLength', nfft);
% figure(4)
% hold on
% for i = 3:14
%     plot(coeffs(:,i))
% end
% hold off


%%
%Pick-off NaN "frameMatchingPhonemes" 4.5
nans = find(isnan(frameMatchingPhonemes)); % find indexes of NaN values
cleanedPhonemeLabels = frameMatchingPhonemes;
cleanedInputMfcc = (mfccCoeff);

cleanedInputMfcc = padarray(cleanedInputMfcc, [0 length(frameMatchingPhonemes) - length(cleanedInputMfcc)], 'post');%DELETE

cleanedPhonemeLabels(nans) = [];
cleanedInputMfcc(:,nans) = [];

%%
%ANN Model DATA SEPERATION 5
tr = [cleanedPhonemeLabels' cleanedInputMfcc'];
n = size(tr, 1);                    % number of samples in the dataset
targets  = tr(:,1);                 % 1st column is |label|
targets(targets == 0) = 10;         % use '10' to present '0'
targetsd = dummyvar(targets);       % convert label into a dummy variable
inputs = tr(:,1:end);               % the rest of columns are predictors

inputs = inputs';                   % transpose input
targets = targets';                 % transpose target
targetsd = targetsd';               % transpose dummy variable

rng(1);                             % for reproducibility
c = cvpartition(n,'Holdout',floor(n/3));   % hold out 1/3 of the dataset

Xtrain = inputs(:, training(c));    % 2/3 of the input for training
Ytrain = targetsd(:, training(c));  % 2/3 of the target for training
Xtest = inputs(:, test(c));         % 1/3 of the input for testing
Ytest = targets(test(c));           % 1/3 of the target for testing
Ytestd = targetsd(:, test(c));      % 1/3 of the dummy variable for testing

%processInputs = {inputs(1,:); inputs(2,:); inputs(3,:); inputs(4,:); inputs(5,:); inputs(6,:); inputs(7,:); 
%                inputs(8,:); inputs(9,:); inputs(10,:); inputs(11,:); inputs(12,:); inputs(13,:);   inputs(14,:);  inputs(15,:); };


%%
clearvars -except inputs targetsd targets processInputs Xtest Ytest
x = inputs(2:14,:);
t = targetsd;


%trainFcn = 'trainscg';  % Scaled conjugate gradient backpropagation.
%trainFcn = 'traingdx';  %90


%net = feedforwardnet(10, trainFcn);
%net.layers{1}.transferFcn = 'poslin';
%net.layers{2}.transferFcn = ''

%%%%%%%%%%%%%%%%%%%%%%%%%
net = layrecnet(1:6,[25]); %patternnet alternative (99.8735 similarity)
%net.trainFcn = 'traincgb';
%net.divideFcn = 'divideblock';
net.layers{2}.transferFcn = 'softmax';
%net.trainFcn = 'trainlm';
%net.trainParam.show = 5;
net.trainParam.max_fail = 80;
net.trainParam.epochs = 3500;


% net.numinputs = 13;
% net.inputConnect = [1 1 1 1 1 1 1 1 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0];
% net = configure(net,x);

view(net)
% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 60/100;
net.divideParam.valRatio = 20/100;
net.divideParam.testRatio = 20/100;


%Annahilate Pre and Post Processing Methods
%Matlab likes to remove constant rows which makes no impact on the
%training(thus reducing the number of weights on the output) we dont want
%this therefore we are removing this behaviour
  net.inputs{1}.processFcns={};
  net.outputs{1}.processFcns = {};
  net.outputs{2}.processFcns={};

  
  
 %p = con2seq(x);
 %pt = con2seq(t);
% Train the Network

 Xgpu = gpuArray(x);
 Tgpu = gpuArray(t);

%LAYRECNET DOES NOT USE GPU ARRAYS
net = configure(net, x, t);
net = train(net, x, t, 'useGPU' , 'yes', 'reduction', 2, 'showResources', 'yes', 'useParallel', 'yes'); %net = train(net, inputs, t, 'useGPU', 'yes');



% Ygpu = net(Xgpu);
% Y = gather(Ygpu);

%[net,trs] = train(net, x, t, 'useParallel', 'yes', 'useGPU', 'yes');


% % Test the Network
% y = net(x);
% e = gsubtract(t,y);
% performance = perform(net,t,y);
% tind = vec2ind(t);
% yind = vec2ind(y);
% %percentErrors = sum(tind ~= yind)/numel(tind)

% View the Network
view(net)

%Write test data
weights_Test_1 = net.IW{1,1}; %net.IW{2,1}, net.b{1}
weights_Test_2 = net.IW;
weights_Test_3 = net.LW{2,1};
net.inputWeights{1,1}.size
wb_2 = getwb(net);
[b,wb_3,wb_4] = separatewb(net, wb_2);
bias1 = net.b{1};
bias2 = net.b{2};

%%
genFunction(net, 'customNetworkFunction_RNN_35K')

% [Xa,Ta] = simpleseries_dataset;
% neta = timedelaynet(1:2,20);
% [Xsa,Xia,Aia,Tsa] = preparets(neta,Xa,Ta);
% neta = train(neta,Xsa,Tsa);
% view(neta)
% Ya = neta(Xsa,Xia,Aia);


%[Xs,Xi,Ai,Ts] = preparets(net,x,t);

Ypred = customNetworkFunction_RNN_35K(Xtest(1:13,:), [], []);     % predicts probability for each label
Ypred(:, 1:5)                      % display the first 5 columns
[~, Ypred] = max(Ypred);            % find the indices of max probabilities
sum(Ytest == Ypred) / length(Ytest) % compare the predicted vs. actual
%%
%Analyze data, compare to sentences
% fileID = fopen('12.txt','r'); %phonesBasic
% ec = textscan(fileID, '%s%s');
% fclose(fileID);
res = net(x(1:13,:), [], []);
[c,cm,ind,pr] = confusion(t, res)


eStart = 1;
eEnd = 800;

%VOA_070202_0630_VOA_Spiker_1-0001 1266, 2.475 8.775
sentence = "polis baghdat'ta düzenlenen digher saldihrihlarda da ondan fazla kishinin öldürüldüghünü bildirdi";

predictions = customNetworkFunction_RNN_35K(x(1:13, eStart:eEnd), [], []); 
[confidence, phonemeIndex] = max(predictions);



B = phonemeIndex(diff([0 phonemeIndex])~=0);
%B = phonemeIndex;

lookupArray = phoneLookup{1}(B);
lookupArray = lookupArray';

strr = "";
for w = 1:length(lookupArray)
    
    chr = lookupArray(w);
    
    if strcmp(chr,'SIL')
        chr = '_';
    elseif strcmp(chr,'ih')
        chr = 'i';
    elseif strcmp(chr,'gh')
        chr = 'g';
    elseif strcmp(chr,'sh')
        chr = 's';
    end
    
    
    
    
    strr = strcat(strr, chr);
end
'done'
%%
% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotconfusion(t,y)
%figure, plotroc(t,y)



%t = phonemesLabelsNumeric;

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
%trainFcn = 'traingdx'; % %90


%net = feedforwardnet(10, trainFcn);
%net.layers{1}.transferFcn = 'poslin';
%net.layers{2}.transferFcn = ''

%%%%%%%%%%%%%%%%%%%%%%%%%
% net = layrecnet(1:2,10, trainFcn); %patternnet alternative (99.8735 similarity)
% net.numinputs = 1;%length(mfccCoeff(:,1));
% %net.numoutputs = length(phoneLookup{1});
% 
% [Xs,Xi,Ai,Ts] = preparets(net,x,t);
% net.layers{2}.transferFcn = 'softmax';
% net.trainParam.max_fail = 12;
% 
% 
% % Setup Division of Data for Training, Validation, Testing
% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio = 15/100;
% net.divideParam.testRatio = 15/100;
% 
% net = layrecnet(1,50);
% net.trainFcn = 'trainbr';
% net.trainParam.show = 5;
% net.trainParam.epochs = 50;
% net.layers{2}.transferFcn = 'softmax';
% view(net);
% net = train(net,x, t);
% 
% y = net(x);
% plot(cell2mat(y))
% 
% %Annahilate Pre and Post Processing Methods
% %Matlab likes to remove constant rows which makes no impact on the
% %training(thus reducing the number of weights on the output) we dont want
% %this therefore we are removing this behaviour
% 
% %net.inputs{1}.processFcns={};
% %net.outputs{1}.processFcns = {};
% %net.outputs{2}.processFcns={};
% 
% 
%  
% % Train the Network
% [net,tr] = train(net,Xs,Ts,Xi,Ai);%train(net,x,t);
% view(net);
% % Test the Network
% y = net(x);
% e = gsubtract(t,y);
% performance = perform(net,t,y)
% tind = vec2ind(t);
% yind = vec2ind(y);
% percentErrors = sum(tind ~= yind)/numel(tind)
% 
% % View the Network
% view(net)
% 
% %Write test data
% weights_Test_1 = net.IW{1,1}; %net.IW{2,1}, net.b{1}
% weights_Test_2 = net.IW;
% weights_Test_3 = net.LW{2,1};
% net.inputWeights{1,1}.size;
% wb_2 = getwb(net);
% [b,wb_3,wb_4] = separatewb(net, wb_2);
% bias1 = net.b{1};
% bias2 = net.b{2};


% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotconfusion(t,y)
%figure, plotroc(t,y)





