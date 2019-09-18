% Simulating protein production and fold change
clear all
close all

% Rep: Reporter
% Flu: Fluorescent reporter

Fracoffs = [0.1 0.05 0.1 0.15 0.2 0.25 0.4];
%Fracoffs = [0.2];

sweepthrough = numel(Fracoffs); %how many values of Fracoff you want to use

minmaxSig = NaN(1,sweepthrough);
minmaxOut = NaN(1,sweepthrough);

for ii = 1:sweepthrough

Phesynthrate = 1; %Reporter synthesis rate. 1/s
Fracoff = Fracoffs(ii);
dil = 1/6000*log(2); % Dilution rate due to volume increase. Assuming doubling time of cell volume = 6000 seconds
Fmat = log(2)/(5.63*60); % sfGFP maturation time. Reference: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=110546

totaltime = 16*60*60; %seconds
dt = 1; %seconds

Rep = zeros(1,totaltime/dt);
Flu = zeros(1,totaltime/dt);
Storesynth = zeros(1,totaltime/dt);

timeaxis = [0:dt:(totaltime-dt)]/3600;

%Phe(1) = 100;

for tt = 2:totaltime/dt
   
    %make synthesis of pheromone change with time as a step function
    
    if mod(tt*dt,100*60) < 33*60 % make pheromone synthesis be a step function so that it is produced for 1/3rd of the cell cycle. 
        Phesynth = Phesynthrate;
    else
        Phesynth = Phesynthrate*Fracoff;
    end
    
    Storesynth(tt) = Phesynth;
    
    Rep(tt) = Rep(tt-1) + dt*(Phesynth - Fmat*Rep(tt-1) - dil*Rep(tt-1));
    Flu(tt) = Flu(tt-1) + dt*(Fmat*Rep(tt-1) - dil*Flu(tt-1));
    
end

%cut down to last two cycles

cutdown = 200*60/dt;

Storesynth = Storesynth((end-cutdown):end)/max(Storesynth((end-cutdown):end));
Rep = Rep((end-cutdown):end)/max(Rep((end-cutdown):end));
Flu = Flu((end-cutdown):end)/max(Flu((end-cutdown):end));
timeaxis = timeaxis((end-cutdown):end) - min(timeaxis((end-cutdown):end));

subplot(3,sweepthrough,0*sweepthrough+1 + ii - 1)
plot(timeaxis,Storesynth)
title([ 'Input signal (Translation rate) 1-TP:' num2str(1-min(Storesynth)/max(Storesynth)) ] );
ylabel('Conc./s')
ylim([0 1.5])
axis square

subplot(3,sweepthrough,1*sweepthrough+1 + ii - 1)
plot(timeaxis,Rep)
title([ 'Reporter 1-TP:' num2str(1-min(Rep)/max(Rep)) ])
ylabel('Conc.')
ylim([0 1.5])
axis square

subplot(3,sweepthrough,2*sweepthrough+1 + ii - 1)
yyaxis left
plot(timeaxis,Storesynth)
ylim([0 1.5])
yyaxis right
plot(timeaxis,Flu)
title([ 'Visible reporter 1-TP:' num2str(1-min(Flu)/max(Flu)) ])
xlabel('Time (hours)')
ylim([0 1.5])
axis square

minmaxSig(ii) = max(Storesynth)/min(Storesynth);
minmaxOut(ii) = 1-min(Flu)/max(Flu);

end

figure
hold on
plot(minmaxSig,minmaxOut,'.','MarkerSize',20)
a = 0.2726;  b = -0.2306 ; %obtained from fitting y = a*(1-exp(b*x)) -- no particular reason to use an exponential other than it fit this data set well
plot([0:0.01:20],a*(1-exp(b*[0:0.01:20])));

plot( [0 log(1-0.2/a)/b] , [0.2 0.2] ,'--k'  )
plot( [log(1-0.2/a)/b log(1-0.2/a)/b] , [0 0.2] ,'--k'  )

xlabel('Input signal FC')
ylabel('Output signal FC')
xlim([0 20])
%ylim([0 0.3])
axis square