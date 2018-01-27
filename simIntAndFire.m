function [simulatedVoltage,meanSpikeRate,instantaneousSpikeRate,spikeVariation,simulatedgSra] = simIntAndFire(neuronParameter,simulatedTime,simulatedInputCurrent)
%% simulatedVoltage,meanSpikeRate,instantaneousSpikeRate,spikeVariation,simulatedgSra] = simIntAndFire(neuronParameter,simulatedTime,simulatedInputCurrent)

% This function simulates the firing of a neuron using the integrate-and-fire model 
% Input:
%       neuronParameter: the structure containing all parameters of the
%                     simulated neuron
%                         neuronParameter.timeConstantMembrane = 0.010;
%                         neuronParameter.VrestingMembrane = -0.065;
%                         neuronParameter.RrestingMembrane = 1e7; % the total membrane resistance
%                         neuronParameter.Vthreshold = -0.050;
%                         neuronParameter.Vreset = -0.065;
%                         neuronParameter.voltageNoise = 0.0000;
%                         neuronParameter.timeConstantSra = 0.10; % time constant for adaptation conductance                         
%                         neuronParameter.Vpotassium = -0.07;     % the simulated resting potential for adaptation
%                         neuronParameter.rm = 0.06;              % the resistance density                        
%                         neuronParameter.spikeDeltaG = 1;        % time step of adaptation update
%       simulatedTime: the timestamp vector of simulated time
%       simulatedInputCurrent: the vector of input current magnitude (must
%                      be the same length as the simulatedTime)
% Output:
%       simulatedVoltage: the resulting membrane voltage
%       meanSpikeRate: the average spike rate of the whole simulation
%       instantaneousSpikeRate: the spike rate at each timestamp (either
%                       calculated or interpolated)
%       spikeVariation: the normalized variation in timing between spikes 
%       simulatedgSra: the conductance of a hypothetical channel
%                       incorporating the refractory pediod

% 1/22/10  dhb  Wrote it as a sub function in integrateAndFireHW.
% 9/9/13   ll   Make it a separate function 

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
simulatedVoltage(1) = neuronParameter.Vreset;
simulatedgSra(1) = 0;
for t = 1:length(simulatedTime)-1
    % Compute instantaneous dV and update simulated voltage.
    dV = (1/neuronParameter.timeConstantMembrane)*(neuronParameter.VrestingMembrane-simulatedVoltage(t) ...
        - neuronParameter.rm*simulatedgSra(t)*(simulatedVoltage(t)-neuronParameter.Vpotassium) ...
        + neuronParameter.RrestingMembrane*simulatedInputCurrent(t))*deltaT;
    if (neuronParameter.voltageNoise > 0)
        dV = dV+neuronParameter.voltageNoise*randn;
    end
    simulatedVoltage(t+1) = simulatedVoltage(t)+dV;
    deltaG = -(1/neuronParameter.timeConstantSra)*simulatedgSra(t)*deltaT;
    simulatedgSra(t+1) = simulatedgSra(t)+deltaG;
    
    
    % Check whether voltage exceeds spike threshold, record a spike, and
    % reset the voltage if so.
    %
    % The spike event bookkeeping here is useful for computing both
    % mean spike rate over the whole train and the time varying
    % spike rate.
    if (simulatedVoltage(t+1) >= neuronParameter.Vthreshold)
        simulatedVoltage(t+1) = neuronParameter.Vreset;
        spikeEvents(t+1) = 1;
        spikeIndex = spikeIndex+1;
        spikeTimes(spikeIndex) = simulatedTime(t+1);
        simulatedgSra(t+1) = simulatedgSra(t+1) + neuronParameter.spikeDeltaG;
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
