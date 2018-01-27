%% Simulate the firing of a simplified neuron using integrate-and-fire
%% model
clear
close all

%% Set up the constants and parameter (in SI unit)
tau_m = 0.01; % the membrane time constant (s)
E_L = -0.065; % the resting potential (V)
Vreset = -0.065; % the reset voltage
Vth = -0.05; % the threshold voltage
V0 = -0.06; % the initial meambrane voltage
Rm = 10^7; % the total membrane resistance (Ohm)
dt = 0.0001; % the time step (s)
tIni = 0; % the initial time
tFin = 0.5; % the final time
numStep = (tFin-tIni)/dt + 1; % number of time steps

timeTrack = tIni:dt:tFin;

%% Simulate the membrane voltage with constant injected current Ie
V = zeros(1,numStep); % initiallize the membrane voltage
V(1) = V0;
Ie = 2e-9; % the injected current (A)

for ii = 2 : numStep
    dV = (1/tau_m) * (E_L - V(ii-1) + Rm * Ie) * dt;
    V(ii) = V(ii-1) + dV;
    if V(ii) >= Vth
        V(ii) = Vreset;
    end
end

% Plot the result
scrPara = get(0,'screensize');
scrWidth = scrPara(3) - scrPara(1);
scrHeight = scrPara(4) - scrPara(2);
h_const = figure('Position',[scrWidth/4 scrHeight/4 scrWidth/2 scrHeight/2]);
plot(timeTrack, V)
title('The simulated membrane voltage of a neuron injected a constant current','FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Membrane voltage (V)','FontSize',15)
set(gca,'ylim',[-0.07 -0.045],'FontSize',15)
V_const = V;

%% Show the firing rates of neuron corresponding to several Ie
Ie = [1:0.1:4] * 10^(-9); % Set several sensible values for Ie
firingRate = zeros(1,length(Ie));
for jj = 1 : length(Ie)
    V = zeros(1,numStep); % initiallize the membrane voltage
    V(1) = V0;
    for ii = 2 : numStep
        dV = (1/tau_m) * (E_L - V(ii-1) + Rm * Ie(jj)) * dt;
        V(ii) = V(ii-1) + dV;
        if V(ii) >= Vth
            firingRate(jj) = firingRate(jj) + 1;
            V(ii) = Vreset;
        end
    end
end
firingRate = firingRate/(tFin - tIni);

% Plot the result
scrPara = get(0,'screensize');
scrWidth = scrPara(3) - scrPara(1);
scrHeight = scrPara(4) - scrPara(2);
figure('Position',[scrWidth/4 scrHeight/4 scrWidth/2 scrHeight/2])
plot(Ie*10^9, firingRate,'r*-')
title('The simulated firing rate with respect to injected current Ie','FontSize',15)
xlabel('Injected current Ie (nA)','FontSize',15)
ylabel('Firing rate (spike/sec)','FontSize',15)
set(gca,'ylim',[min(firingRate)-10 max(firingRate)+10],'FontSize',15)

%% Simulate the dynamics of membrane voltage corresponding to fluctuating current Ie
% Set two waveforms for fluctutating current 
% Sinusoidal wave with period T1
T1 = 0.05;
Ie1 = sin(2*pi*timeTrack/T1) * 4e-9; 

% Sawtooth wave with period T2
T2 = 0.07;
Ie2 = sinc(2*pi*(timeTrack-max(timeTrack)/2)/T2) * 4e-9;

% Combine the two waveforms
Ie = [Ie1;Ie2];

% Initiallize the membrane voltage
V = zeros(2,numStep); 
V(:,1) = [V0,V0];

% Simulate the membrane voltage with fluctuating injected current Ie
for jj = 1 : 2
    for ii = 2 : numStep
        dV = (1/tau_m) * (E_L - V(jj,ii-1) + Rm * Ie(jj,ii)) * dt;
        V(jj,ii) = V(jj,ii-1) + dV;
        if V(jj,ii) >= Vth
            V(jj,ii) = Vreset;
        end
    end
end

% Plot the result
scrPara = get(0,'screensize');
scrWidth = scrPara(3) - scrPara(1);
scrHeight = scrPara(4) - scrPara(2);
figure('Position',[scrWidth/7 scrHeight/7 5*scrWidth/7 5*scrHeight/7])
title('The simulated membrane voltage of a neuron','FontSize',15)

subplot(2,2,1)
plot(timeTrack, V(1,:))
title('The simulated membrane voltage of a neuron','FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Membrane voltage (V)','FontSize',15)
set(gca,'ylim',[min(V(1,:))-0.01  max(V(1,:)) + 0.01],'FontSize',15)

subplot(2,2,2)
plot(timeTrack, V(2,:))
title('The simulated membrane voltage of a neuron','FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Membrane voltage (V)','FontSize',15)
set(gca,'ylim',[min(V(2,:))-0.01  max(V(2,:))+0.01],'FontSize',15)

subplot(2,2,3)
plot(timeTrack, Ie(1,:)*10^9)
title('The simulated injected current of a neuron','FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Injected current (nA)','FontSize',15)
set(gca,'ylim',[min(Ie(1,:)*10^9)-1  max(Ie(1,:)*10^9)+1],'FontSize',15)

subplot(2,2,4)
plot(timeTrack, Ie(2,:)*10^9)
title('The simulated injected current of a neuron','FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Injected current (nA)','FontSize',15)
set(gca,'ylim',[min(Ie(2,:)*10^9)-1  max(Ie(2,:)*10^9)+1],'FontSize',15)
    
%% Simulate the membrane voltage with constant injected current Ie and Gaussian 
%% random noise induced at the membrane voltage
V = zeros(1,numStep); % initiallize the membrane voltage
V(1) = V0;
Ie = 2e-9; % the injected current (A)

for ii = 2 : numStep
    dV = (1/tau_m) * (E_L - V(ii-1) + Rm * Ie) * dt;
    V(ii) = V(ii-1) + dV;
    if V(ii) >= Vth
        V(ii) = Vreset;
    end
    % Add some Gaussian noise of mean 0 and std 0.001
    V(ii) = V(ii) + randn * 0.001;
end

% Plot the result
figure(h_const)
hold on
plot(timeTrack, V,'r')
legend('Noise-free membrane voltage','Noisy membrane voltage')

%% Simulate the membrane voltage with constant Ie and incorporated
%% refractory period
% Set up the parameters
tau_sra = 0.1; % time constant for adaptation conductance
rmTimeGsra = 0; % initial value of adaptation conductance
rmTimeDelgsra = 0.06;
Ek = -0.07;
V = zeros(1,numStep); % initiallize the membrane voltage
V(1) = V0;
Ie = 2e-9; % the injected current (A)

for ii = 2 : numStep
    dV = (1/tau_m) * (E_L - V(ii-1) + Rm * Ie - rmTimeGsra * (V(ii-1)-Ek)) * dt;
    V(ii) = V(ii-1) + dV;
    if V(ii) >= Vth
        V(ii) = Vreset;
        rmTimeGsra = rmTimeGsra + rmTimeDelgsra; 
    elseif rmTimeGsra > 0
        dRmTimeGsra = -rmTimeGsra * dt / tau_sra;
        rmTimeGsra = rmTimeGsra + dRmTimeGsra;
    end    
end

% Plot the result
scrPara = get(0,'screensize');
scrWidth = scrPara(3) - scrPara(1);
scrHeight = scrPara(4) - scrPara(2);
figure('Position',[scrWidth/4 scrHeight/4 scrWidth/2 scrHeight/2])
plot(timeTrack, V,'r')
hold on
plot(timeTrack,V_const,'b')
title('Membrane voltage with and without refractory parameter','FontSize',15)
xlabel('Time (s)','FontSize',15)
ylabel('Membrane voltage (V)','FontSize',15)
set(gca,'ylim',[-0.07 -0.045],'FontSize',15)
legend('With refractory parameter','Without refractory parameter')
    