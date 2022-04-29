clear;close all
rng default% reproducibility

%% set up simulation
fs = 30;% sampling rate (Hz)
spdMu = 5;% mm/s
spdSigma = 0.5;
curvMu = 60;% degrees/s
curvSigma = 6;% degrees/s

t2Cons = divisors(fs); % trajectory durations used simulations (in frames)
tHalfSecondDur = find(t2Cons==(fs/2)); % get index of duration=0.5s for plotting purposes

nAgents = 1000; % number of agents
simDuration = 5;% seconds

%% simulate
for i = 1:numel(t2Cons)
    % calculate number of trajectory samples there are given total time and
    % state persistence
    t = t2Cons(i);
    nSamples = simDuration.*fs./t;
    
    % initialize
    Ang = zeros(1000,t*nSamples);
    xModel = zeros(1000,t*nSamples);
    yModel = zeros(1000,t*nSamples);
    
    % Loop through each trajectory
    for tt = 1:nSamples
        initNdx = (tt-1).*t+1;
        
        % sample speed and curvature from distributions
        spd = normrnd(spdMu,spdSigma,1000,1)./fs;
        curv = normrnd(curvMu,curvSigma,1000,1)./fs;
        spd = repmat(spd,1,t);
        curv = repmat(curv,1,t+1);
        
        % update positions of each agent for this trajectory
        for j = 1:nAgents
            % update the orientation of movement during this trajectory
            Ang(j,initNdx:initNdx+t) = Ang(j,initNdx)+[0 (cumsum((curv(j,1:end-1)+curv(j,2:end))./2))];
            
            % update the positions of movement during this trajectory
            xModel(j,initNdx:initNdx+t) = xModel(j,initNdx)+[0 cumsum(spd(j,1:end).*cosd(Ang(j,initNdx:initNdx+t-1)))];
            yModel(j,initNdx:initNdx+t) = yModel(j,initNdx)+[0 cumsum(spd(j,1:end).*sind(Ang(j,initNdx:initNdx+t-1)))];
        end
    end
    
    % fit the final positions to a bivariate gaussian
    GMModel = fitgmdist([xModel(:,end), yModel(:,end)],1);
    
    % get the generalized variance
    genVar(i) = det(GMModel.Sigma);%prod(eig(GMModel.Sigma));
    
    % keep the data from the trial used for the figures
    if i == tHalfSecondDur
        xModel_plotting = xModel;
        yModel_plotting = yModel;
    end
    
end

%% plotting

% plot the generalized variance as a function of trajectory duration
figure;plot(t2Cons./fs,genVar,'Linewidth',1);
ylim([0 4])
xlabel('State duration (seconds)');ylabel('generalized variance (mm^2)')
title('after 5 seconds')

% plot the spread of flies at the half second mark and the one second mark
figure;set(gcf,'Position',[2 42 838 924])
subplot(3,1,1);
scatter(xModel_plotting(:,t2Cons(tHalfSecondDur)),...
    yModel_plotting(:,t2Cons(tHalfSecondDur)))
axis([0 6 0 4]);axis square
xlabel('x-position (mm)');ylabel('y-position (mm)')
title('Position at 0.5 seconds')
subplot(3,1,2);scatter(xModel_plotting(:,fs),yModel_plotting(:,fs))
axis([0 6 0 4]);axis square
xlabel('x-position (mm)');ylabel('y-position (mm)')
title('Position at 1 second')
subplot(3,1,3);
scatter(xModel_plotting(:,t2Cons(tHalfSecondDur)),...
    yModel_plotting(:,t2Cons(tHalfSecondDur)));hold on;
scatter(xModel_plotting(:,fs),yModel_plotting(:,fs))
axis([0 6 0 4]);axis square
xlabel('x-position (mm)');ylabel('y-position (mm)')
title('Position at 0.5 and 1 second')

% plot the speed and curvature distributions
figure;set(gcf,'Position',[2 42 838 924])
subplot(2,1,1);
x = [0:0.1:10];
y = normpdf(x,spdMu,spdSigma);
plot(x,y)
xlabel('Speed');ylabel('PDF')
ylim([0 1])
subplot(2,1,2);
x = [30:0.5:90];
y = normpdf(x,curvMu,curvSigma);
plot(x,y)
ylim([0 0.1])
xlabel('Curvature');ylabel('PDF')

% print to ps file
for f = 1:get(gcf,'number')
    figure(f);
    print('-painters','-dpsc2',['SensorimotorNoiseSimulations.ps'],'-loose','-append');
end

% covert from ps to pdf (need ghostscript for this)
ps2pdf('psfile', 'SensorimotorNoiseSimulations.ps', 'pdffile', ...
    'SensorimotorNoiseSimulations.pdf', 'gspapersize', 'letter',...
    'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
    'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
    'gslibpath','C:\Program Files\gs\gs9.50\lib');

