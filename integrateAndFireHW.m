function integrateAndFireHW
% integrateAndFireHW
%
% Solution to the integrate and fire model homeowork for
% Vijay's comp neuro course.
%
% 1/22/10  dhb  Wrote it.

%% Clear and close
clear; close all;

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
nTimeSteps = 36000;
secondsPerStep = 0.5e-4;

%% Initialize timebase
simulatedTime = (1:nTimeSteps)*secondsPerStep;

%% Simulate for one input current and graph.
%
% If the voltage noise is set to zero, then the spike train is
% nice and regular, the mean of the instaneous spike rate matches
% the separately computed mean spike rate with zero standard deviation,
% and the returned spike variation is 0.
%
% Set neuron.spikeDeltaG = 0 to get rid of adaptation.
inputCurrent0 = 1.6e-9;
simulatedInputCurrent = inputCurrent0*ones(nTimeSteps,1);
[simulatedVoltage,spikeRate,instantaneosSpikeRate,spikeVariation,simulatedgSra] = simIntAndFire(neuron,simulatedTime,simulatedInputCurrent);
fprintf('Input current %0.1f nanoamps, spike rate %0.1f spikes/sec\n',inputCurrent0*1e9,spikeRate);
fprintf('Mean of instantaneous spike rate %0.1f +/- %0.1f spikes/sec\n',...
    mean(instantaneosSpikeRate(~isnan(instantaneosSpikeRate))),std(instantaneosSpikeRate(~isnan(instantaneosSpikeRate))));
fprintf('Returned spike time variation is %0.1f\n',spikeVariation);
figure; clf;
subplot(2,1,1); hold on
plot(simulatedTime*1000,simulatedVoltage*1000,'r');
xlabel('Time (mSec)');
ylabel('Voltage (mV)');
ylim([-75 -45]);
xlim([0 max(simulatedTime*1000)]);
subplot(2,1,2); hold on
plot(simulatedTime*1000,simulatedgSra,'r');
xlabel('Time (mSec)');
ylabel('gSRA)');
xlim([0 max(simulatedTime*1000)]);

%% Simulate for a range of input currents and plot spike rate versus current
nInputCurrents = 100;
lowInputCurrent = 1e-9;
highInputCurrent = 7e-9;
theInputCurrents = linspace(lowInputCurrent,highInputCurrent,nInputCurrents);
rawSimulatedInputCurrent = ones(nTimeSteps,1);
for i = 1:nInputCurrents
    [nil,spikeRates(i)] = simIntAndFire(neuron,simulatedTime,theInputCurrents(i)*rawSimulatedInputCurrent); %#ok<*AGROW>
end
figure; clf; hold on
plot(theInputCurrents*1e9,spikeRates,'ro','MarkerSize',2','MarkerFaceColor','r');
plot(theInputCurrents*1e9,spikeRates,'r');
xlabel('Input Current (nanoamps)');
ylabel('Spikes/Sec');

%% Simulate spike variation as a function of spike noise
inputCurrent0 = 2e-9;
simulatedInputCurrent = inputCurrent0*ones(nTimeSteps,1);
nNoises = 100;
lowNoise = 0.000001;
highNoise = 0.0005;
theNoises = logspace(log10(lowNoise),log10(highNoise),nNoises);
for i = 1:nInputCurrents
    neuron.voltageNoise = theNoises(i);
    [nil,spikeRates(i),nil,spikeVariation(i)] = simIntAndFire(neuron,simulatedTime,simulatedInputCurrent); %#ok<*AGROW>
end
figure; clf;
subplot(1,2,1); hold on
plot(theNoises*1000,spikeRates,'ro','MarkerSize',2','MarkerFaceColor','r');
plot(theNoises*1000,spikeRates,'r');
xlabel('Voltage Noise Stdev (mV)');
ylabel('Spikes/Sec');
xlim([0 highNoise*1000]);
ylim([0 max(spikeRates)]);
subplot(1,2,2); hold on
plot(theNoises*1000,spikeVariation,'ro','MarkerSize',2','MarkerFaceColor','r');
plot(theNoises*1000,spikeVariation,'r');
xlabel('Voltage Noise Stdev (mV)');
ylabel('Interspike Variation');
xlim([0 highNoise*1000]);
ylim([0 max(spikeVariation)]);

%% Simulate for a time varying input currents and plot
% time varying spike rate as well input and voltage
% versus time.
timeVaryingType = 'sinusoid';
lowInputCurrent = 2e-9;
highInputCurrent = 5e-9;
frequency = 1;
neuron.voltageNoise = 0;
switch(timeVaryingType)
    case 'sinusoid'
        simulatedInputCurrent = lowInputCurrent+(highInputCurrent-lowInputCurrent)*(sin(2*pi*frequency*simulatedTime)+1)/2;
    case 'sawtooth'
        simulatedInputCurrent = lowInputCurrent+(highInputCurrent-lowInputCurrent)*(sawtooth(2*pi*frequency*simulatedTime)+1)/2;
    otherwise
        error('Oops. Unknown time varying signal speified\n');
end
[simulatedVoltage,meanSpikeRate,instantaneousSpikeRate] = simIntAndFire(neuron,simulatedTime,simulatedInputCurrent);
figure; clf;
subplot(3,1,1); hold on
plot(simulatedTime*1000,simulatedInputCurrent*1e9,'r');
xlabel('Time (mSec)');
ylabel('Input Current (nanoamps)');
xlim([0 max(simulatedTime*1000)]);
ylim([0 max(simulatedInputCurrent*1e9)]);
subplot(3,1,2); hold on
plot(simulatedTime*1000,simulatedVoltage*1000,'r');
xlabel('Time (mSec)');
ylabel('Voltage (mV)');
xlim([0 max(simulatedTime*1000)]);
ylim([-75 -45]);
subplot(3,1,3); hold on
plot(simulatedTime*1000,instantaneousSpikeRate,'r');
xlabel('Time (mSec)');
ylabel('Spikes/Sec');
xlim([0 max(simulatedTime*1000)]);
ylim([0 max(instantaneousSpikeRate)]);


end


%% [simulatedVoltage,spikeRate,spikeTimes]= simIntAndFire(neuron,simulatedTime,simulatedInputCurrent,voltageNoise)
function [simulatedVoltage,meanSpikeRate,instantaneousSpikeRate,spikeVariation,simulatedgSra] = simIntAndFire(neuron,simulatedTime,simulatedInputCurrent)
%
% This function actually simulates out the model.  

% Check inputs
if (length(simulatedTime) ~= length(simulatedInputCurrent))
    error('Oops.  Length of time and input current vectors must match\n');
end

% Simulate the integrate and fire process

spikeEvents = NaN*zeros(size(simulatedTime));
deltaT = simulatedTime(2)-simulatedTime(1);

spikeIndex = 0;
spikeTimes = [];

% Loop through time to simulate the differential equation.
simulatedVoltage = NaN*zeros(size(simulatedTime));
simulatedgSra = NaN*zeros(size(simulatedTime));
simulatedVoltage(1) = neuron.Vreset;
simulatedgSra(1) = 0;
for t = 1:length(simulatedTime)-1
    % Compute instantaneous dV and update simulated voltage.
    dV = (1/neuron.timeConstantMembrane)*(neuron.VrestingMembrane-simulatedVoltage(t) ...
        - neuron.rm*simulatedgSra(t)*(simulatedVoltage(t)-neuron.Vpotassium) ...
        + neuron.RrestingMembrane*simulatedInputCurrent(t))*deltaT;
    if (neuron.voltageNoise > 0)
        dV = dV+neuron.voltageNoise*randn;
    end
    simulatedVoltage(t+1) = simulatedVoltage(t)+dV;
    deltaG = -(1/neuron.timeConstantSra)*simulatedgSra(t)*deltaT;
    simulatedgSra(t+1) = simulatedgSra(t)+deltaG;
    
    
    % Check whether voltage exceeds spike threshold, record a spike, and
    % reset the voltage if so.
    %
    % The spike event bookkeeping here is useful for computing both
    % mean spike rate over the whole train and the time varying
    % spike rate.
    if (simulatedVoltage(t+1) >= neuron.Vthreshold)
        simulatedVoltage(t+1) = neuron.Vreset;
        spikeEvents(t+1) = 1;
        spikeIndex = spikeIndex+1;
        spikeTimes(spikeIndex) = simulatedTime(t+1);
        simulatedgSra(t+1) = simulatedgSra(t+1) + neuron.spikeDeltaG;
    else
        spikeEvents(t+1) = 0;
    end
end

% Check output
if (any(isnan(simulatedVoltage)) || any (isinf(simulatedVoltage)))
    error('Oops.  Simulation went off the rails\n');
end

% Compute mean spike rate over the whole spike train
spikeIndices = find(spikeEvents == 1);
nSpikes = length(spikeIndices);
if (nSpikes > length(simulatedVoltage)/2)
    fprintf('Warning.  Spike rate is high relative to size of simulated time step\n');
end
if (nSpikes <= 1)
    meanSpikeRate = 0;
else
    measurementInterval = simulatedTime(spikeIndices(end))-simulatedTime(spikeIndices(1));
    meanSpikeRate = (nSpikes-1)/measurementInterval;
end

% Compute spike rate at time of each spike and then inteprolate over the
% time range where there is data for it to make sense to do so.
%
% Also compute mean and variance of interspike intervals.
if (length(spikeTimes) > 3)
    for i = 1:length(spikeTimes)-2;
        rawInstanteousSpikeTime(i) = spikeTimes(i+1);
        rawInstanteousSpikeRate(i) = 1/((spikeTimes(i+2)-spikeTimes(i))/2);
    end
    instantaneousSpikeRate = interp1(rawInstanteousSpikeTime,rawInstanteousSpikeRate,simulatedTime,'linear',NaN);
    spikeIntervals = diff(spikeTimes);
    spikeVariation = std(spikeIntervals)/mean(spikeIntervals);    
else
    instantaneousSpikeRate = [];
    spikeVariation = [];
end


end



