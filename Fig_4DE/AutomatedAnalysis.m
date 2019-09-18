redofits = 1; %decide if you want to do all the fits again

if redofits == 1
clear all
close all

%analyze profiles and fit them to exponential curves

load('Rec_Profiles.mat');

diffusion_constants = [];

pixtomic = 0.0641; % 1 pixel=0.0641 microns

finaltimepoint = 41;

storeallydata = []; %just so you can compare these later

background = 110;

for i = 1:28

frap = allprofiles{2,i}; %get avg intensity in frapped region
wholecell = allprofiles{3,i}; %get intensity of whole cell to correct for bleaching
radius = allprofiles{4,i}*pixtomic/2; %divide length of frapped region by 2 to get radius
time = (0:size(frap,1)-1)*15;

%normalize fluoresence
fluor = (frap-background)./(wholecell-background);
%fluor = frap;
%% carry out fit
[xData, yData] = prepareCurveData( time, fluor );

% Set up fittype and options.
%ft = fittype( 'c-a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a*(1-exp(-b*x))+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.Robust = 'LAR';
%opts.StartPoint = [0.777273516432061 0.639333731520255 0.827697261051344];
opts.StartPoint = [1 1 1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

halflife = log(2)/fitresult.b;

D = 0.22*radius^2/halflife;

diffusion_constants = [diffusion_constants D];

%% Generate plots

subplot(4,7,i)
hold on% Plot fit with data.
h = plot( fitresult, xData, yData );
%legend( h, 'Data', 'Fit', 'Location', 'SouthEast' );
legend off
% Label axes
xlabel Time(s); ylabel('Intensity');
xlim([0 630]); ylim([0 1]);
axis square


storeallydata = [storeallydata , yData] ;

end

end

figure

scatter(1+0.1*rand(1,5),diffusion_constants(1:5),'fill');
hold on
scatter(1+0.3*rand(1,2),diffusion_constants(6:7),'fill');
scatter(1+0.3*rand(1,7),diffusion_constants(8:14),'fill');
scatter(1+0.3*rand(1,3),diffusion_constants(15:17),'fill');
scatter(2+0.3*rand(1,5),diffusion_constants(18:22),'fill');
scatter(2+0.3*rand(1,6),diffusion_constants(23:28),'fill');

plot( [0.75 1.25] , [ mean(diffusion_constants(1:17)) mean(diffusion_constants(1:17)) ] );
plot( [1.75 2.25] , [ mean(diffusion_constants(18:28)) mean(diffusion_constants(18:28)) ] );

xlim([0 3])

mean(diffusion_constants(1:17))

mean(diffusion_constants(18:28))


%generate example graphs
figure
for i = [2,19]

hold on

yData=storeallydata(:,i);
% Set up fittype and options.
%ft = fittype( 'c-a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a*(1-exp(-b*x))+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf];
opts.Robust = 'LAR';
%opts.StartPoint = [0.777273516432061 0.639333731520255 0.827697261051344];
opts.StartPoint = [1 1 1];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
plot( fitresult, xData, yData );
%plot(xData,yData,'.');

end

