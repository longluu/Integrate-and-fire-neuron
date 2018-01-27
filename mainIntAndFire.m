%% Run the simulation of neuron's firing from function simIntAndFire which employs the integrate-and-fire model
%% Clear and close
clear; close all;
clc

%% Put analysis routines on path
PupilAnalysisPath = '/Users/Shared/Matlab/Experiments/OneLight/OLPupilDiameter/analysis/OLSequentialTrialAnalysisFunctions';
if (~any(findstr(path,'OLSequentialTrialAnalysisFunctions')))
    addpath(genpath(PupilAnalysisPath),'-end');
end

%% Parameters
neuron.timeConstantMembrane = 0.010;
neuron.VrestingMembrane = -0.065;
neuron.RrestingMembrane = 1e7;
neuron.Vthreshold = -0.050;
neuron.Vreset = -0.065;
neuron.voltageNoise = 0.0000;
neuron.timeConstantSra = 0.10;
neuron.Vpotassium = -0.07;
neuron.rm = 0.06;
neuron.spikeDeltaG = 1;
neuron.voltageNoise = 0;
lowInputCurrent = 2e-9;
highInputCurrent = 5e-9;
frequency = [0.01 0.05 0.1 0.5 1 2];
nSimulateCycles = 2;
simulateTimestep = 1e-4;
simulateEndtime = nSimulateCycles * (1./frequency);

%% Initialize timebase
simulateTime = cell(1,length(frequency));
for ii = 1 : length(frequency)
    simulateTime{ii} = 0:simulateTimestep:simulateEndtime(ii);
end

%% Estimate the temporal transfer function
ampFit = NaN(1,length(frequency));
phaseFit = NaN(1,length(frequency));
spikeRateFit = cell(1,length(frequency));
maxSpike = 0;
simulatedInputCurrent = cell(1,length(frequency));
instantaneousSpikeRate = cell(1,length(frequency));
simulatedTimeFit = cell(1,length(frequency));
simulatedVoltage = cell(1,length(frequency));
meanSpikeRate = NaN(1,length(frequency));
maxRate = 0;
for ii = 1 : length(frequency)
    simulatedInputCurrent{ii} = lowInputCurrent+(highInputCurrent-lowInputCurrent)*(sin(2*pi*frequency(ii)*simulateTime{ii})+1)/2;
    [simulatedVoltage{ii},meanSpikeRate(ii),instantaneousSpikeRate{ii}] = simIntAndFire(neuron,simulateTime{ii},simulatedInputCurrent{ii});
    matchParameter.sgolay_span = 20; % support of sgolay impulse response
    matchParameter.sgolay_polynomial = 7; % degree of sgolay polynomial
    matchParameter.sampling_frequency = 20; 
    simulatedTimeMs = simulateTime{ii} * 1000; % convert to ms for compatibility with SGgolaySmoothTest function
    [ampFit(ii), phaseFit(ii), spikeRateFit{ii}, simulatedTimeFit{ii}] = SpectralResponse(simulatedTimeMs, instantaneousSpikeRate{ii}, frequency(ii), matchParameter);
    maxRate = max([maxRate spikeRateFit{ii}]);
end

%% Plot the input current, spike rate and fitted spike rate
for ii = 1 : length(frequency)
    scrPara = get(0,'screensize');
    scrWidth = scrPara(3) - scrPara(1);
    scrHeight = scrPara(4) - scrPara(2);
    figure('Position',[scrWidth/3 scrHeight/3 scrWidth/3 scrHeight/3])
    fontSize = 15;
    
    subplot(2,1,1); hold on
    set(gca,'FontSize',15)
    plot(simulateTime{ii},simulatedInputCurrent{ii}*1e9,'r');
    plot(simulateTime{ii}, mean([lowInputCurrent,highInputCurrent])*1e9 * ones(1,length(simulateTime{ii})), 'b-') 
    [maxInput, maxInd] = max(simulatedInputCurrent{ii}(1:round(length(simulatedInputCurrent{ii})/2))*1e9);
    plot([simulateTime{ii}(maxInd) simulateTime{ii}(maxInd)], [0 maxInput], 'g');
    xlabel('Time (Sec)');
    ylabel('Input Current (nanoamps)');
    xlim([0 max(simulateTime{ii})]);
    ylim([0 max(simulatedInputCurrent{ii}*1e9)]);
    title(['Input frequency ' num2str(frequency(ii))]);
    
    subplot(2,1,2); hold on
    set(gca,'FontSize',15)
    plot(simulateTime{ii},instantaneousSpikeRate{ii},'r');
    xlabel('Time (Sec)');
    ylabel('Spike rate (spikes/sec)');
    xlim([0 max(simulateTime{ii})]);
    ylim([0 maxRate + 50]);
    plot(simulatedTimeFit{ii},spikeRateFit{ii},'b');
    [maxInput, maxInd] = max(spikeRateFit{ii}(1:round(length(spikeRateFit{ii})/2)));
    plot([simulatedTimeFit{ii}(maxInd) simulatedTimeFit{ii}(maxInd)], [0 maxInput], 'g');        
    legend('Simulated spike rate','Fitted spike rate')
    savefig(['Frequency_' num2str(ii)],gcf,'pdf')
%     subplot(3,1,3);
%     set(gca,'FontSize',15)
%     plot(simulateTime{ii},simulatedVoltage{ii},'r');
%     xlabel('Time (Sec)');
%     ylabel('Simulated volage (V)');
%     xlim([0 max(simulateTime{ii})]);
end

%% Plot the transfer function
scrPara = get(0,'screensize');
scrWidth = scrPara(3) - scrPara(1);
scrHeight = scrPara(4) - scrPara(2);
figure('Position',[scrWidth/3 scrHeight/3 scrWidth/3 scrHeight/3])
set(gca,'FontSize',15)
sizeMarker = 10;
semilogx(frequency, ampFit, 'bo-','MarkerSize',sizeMarker)
title('The amplitude of transfer function')
xlabel('Frequency (logHz)')
ylabel('Amplitude')
axis([frequency(1)-0.5 frequency(end)+0.5 0 max(ampFit)+1])
savefig('MTF_amplitude',gcf,'pdf')
    
scrPara = get(0,'screensize');
scrWidth = scrPara(3) - scrPara(1);
scrHeight = scrPara(4) - scrPara(2);
figure('Position',[scrWidth/3 scrHeight/3 scrWidth/3 scrHeight/3])
set(gca,'FontSize',15)
semilogx(frequency, phaseFit, 'bo-','MarkerSize',sizeMarker)
title('The phase of transfer function')
xlabel('Frequency (logHz)')
ylabel('Phase (radian)')
axis([frequency(1)-0.5 frequency(end)+0.5 min(phaseFit)-0.5 max(phaseFit)+0.5])
savefig('MTF_phase',gcf,'pdf')
