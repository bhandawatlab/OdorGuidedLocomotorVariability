clear;close all
rng default% reproducibility

%% set up simulation
nAgent = 5000; % number of agents/simulation
nTrials = 21; % number of simulations
fs = 30; % sampling rate (frames/s)
simulation_duration = 120;% seconds
tt = simulation_duration*fs;

% statistics of environment
plumeRate = 1./(fs*2);% hz
plumeSD = 1./(fs*10);% hz
dutyCycle = 20;
t = 0:1:tt;
grids = [-600:10:600];

% set up agent parameters
onSpd = 10;% mm/s
% each simulation uses a different off speed 
offSpd = onSpd-([1:nTrials]-1)./2;% mm/s
spdSigma = 0;
curvMu = 60;%degrees/s
curvSigma = 0;%
corr = 0.5;% walk correlation with previous direction. 0.5 means uncorrelated
dur = fs./2;% duration in steps

% get index where on speed is half of off speed for plotting purposes
tHalfSpeed = find(offSpd==(onSpd/2)); 


%% simulate
stimMeanAllTrial = zeros(nAgent,nTrials);
stimVarAllTrial = zeros(nAgent,nTrials);
for trial = 1:nTrials
    
    % set up stimulus environment
    [XX,YY] = meshgrid(grids,grids);
    y = (square(2*pi*(plumeRate+(randn(numel(XX),1).*plumeSD))*(0:1:tt),dutyCycle)+1)./2;
    
    % prealocate
    Ang = zeros(nAgent,tt);
    xModel = zeros(nAgent,tt);
    yModel = zeros(nAgent,tt);
    dirPrev = ones(nAgent,1);
    spdOffSampled = normrnd(offSpd(trial),spdSigma,nAgent,1)./fs;%10.*speedFact(trial)
    curvSampled = normrnd(curvMu,curvSigma,nAgent,1)./fs;
    curvAll = zeros(nAgent,tt).*curvSampled;
    spdAll = zeros(nAgent,tt).*spdOffSampled;
    
    for t = 2:tt
        initNdx = (tt-1).*t+1;
        if mod(t,dur) == 1
            % sample from distribution
            spdOnSampled = normrnd(onSpd,spdSigma,nAgent,1)./fs;
            spdOffSampled = normrnd(offSpd(trial),spdSigma,nAgent,1)./fs;%10.*speedFact(trial)
            curvSampled = normrnd(curvMu,curvSigma,nAgent,1)./fs;
            
            % calculate direction
            dir = dirPrev.*(2*(randn(nAgent,1)<=corr)-1);
            dirPrev = dir;
            curvSampled = curvSampled.*dir;
            
            % set on new speed as on or off based on last sensory experience
            [~,~,~,binX,binY] = histcounts2(xModel(:,t-1),yModel(:,t-1),grids,grids);
            ind = sub2ind(size(XX),binX+1,binY+1);
            spdOffSampled(y(ind,t)==1) = spdOnSampled(y(ind,t)==1);
        end
        
        % update positions of each agent
        curvAll(:,t) = curvSampled;
        spdAll(:,t) = spdOffSampled;
        Ang(:,t) = Ang(:,t-1)+(curvAll(:,t-1)+curvAll(:,t))./2;
        xModel(:,t) = xModel(:,t-1)+spdAll(:,t).*cosd(Ang(:,t));
        yModel(:,t) = yModel(:,t-1)+spdAll(:,t).*sind(Ang(:,t));
    end

    % Calculate the stimulation experienced by each agent based on position
    [~,~,~,binX,binY] = histcounts2(xModel,yModel,grids,grids);
    ind = sub2ind(size(XX),binX+1,binY+1);
    stim = zeros(size(ind));
    for t = 1:tt
        stim(:,t) = y(ind(:,t),t);
    end
    
    % keep the data from the trial used for the figures
    if trial==tHalfSpeed
        xModel_plotting = xModel;
        yModel_plotting = yModel;
        spdAll_plotting = spdAll;
        stim_plotting = stim;
    end
    
    % calculate the mean and variance of the stimulation experienced across
    % time for each agent
    stimMeanAllTrial(:,trial) = mean(stim,2);
    stimVarAllTrial(:,trial) = var(stim,[],2);
end

%% plotting

% plot the histogram of the stimulus mean and variance experienced across
% two batches of flies (agents) 1.) when no change in speed, and 2.) when
% speed is 0
figure;set(gcf,'Position',[2 42 838 924]);
subplot(2,2,1);
histogram(stimMeanAllTrial(:,1),[0:0.005:0.4],'Normalization','Probability');hold on;
histogram(stimMeanAllTrial(:,end),[0:0.005:0.4],'Normalization','Probability')
legend({'faster movement','slower movement'})
xlabel('Stimulus average over time');ylabel('Probability')
title('Mean');
subplot(2,2,2);
histogram(stimVarAllTrial(:,1),[0:0.005:0.3],'Normalization','Probability');hold on;
histogram(stimVarAllTrial(:,end),[0:0.005:0.3],'Normalization','Probability')
legend({'faster movement','slower movement'})
xlabel('Stimulus variance over time');ylabel('Probability')
title('Variance')

% plot the mean of the experienced stimulus mean and mean of the stimulus
% variance
subplot(2,2,3);
plot(10-[0:(nTrials-1)]./2,mean(stimMeanAllTrial));
xlabel('speed at off (mm/s)');ylabel('Mean stimulus mean')
%ylim([0.2 0.23])
subplot(2,2,4);
plot(10-[0:(nTrials-1)]./2,mean(stimVarAllTrial));
xlabel('speed at off (mm/s)');ylabel('Mean stimulus variance')
%ylim([0.155 0.175])

% set up sample time frames to plot
tLim = 600:900;flySamp = 1;
tt = (0:range(tLim))./fs;

% plot the x,y paths
figure;set(gcf,'Position',[2 42 838 924])
subplot(8,1,[1:4]);
plot(xModel_plotting(flySamp,tLim),yModel_plotting(flySamp,tLim));hold on;
plot(xModel_plotting(flySamp,tLim(1)),yModel_plotting(flySamp,tLim(1)),'r*');
xticks(grids)
yticks(grids)
grid on
axis equal
xlabel('x-position (mm)');ylabel('y-position (mm)')

% plot the stimulus in two sample grids over the sample time frame
subplot(8,1,5);
plot(tt,y(1,tLim))
xlim([0 max(tt)]);ylim([0 1])
xlabel('time (s)');title('Stimulus at sample square 1')
subplot(8,1,6);
plot(tt,y(2,tLim))
xlim([0 max(tt)]);ylim([0 1])
xlabel('time (s)');title('Stimulus at sample square 1')

% plot the stimulus experienced by the agent over the sample time frame
subplot(8,1,7);
plot(tt,stim_plotting(flySamp,tLim))
xlim([0 max(tt)]);ylim([0 1])
xlabel('time (s)');title('Stimulus experienced by fly')
% plot the speed profile by the agent over the sample time frame
subplot(8,1,8);
plot(tt,spdAll_plotting(flySamp,tLim).*fs)
xlim([0 max(tt)]);ylim([0 onSpd])
xlabel('time (s)');title('Speed performed by fly')

% print to ps file
for f = 1:get(gcf,'number')
    figure(f);
    print('-painters','-dpsc2',['EnvironmentMovementSimulations.ps'],'-loose','-append');
end

% covert from ps to pdf (need ghostscript for this)
ps2pdf('psfile', 'EnvironmentMovementSimulations.ps', 'pdffile', ...
    'EnvironmentMovementSimulations.pdf', 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');






