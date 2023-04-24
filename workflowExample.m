
clear; close all; clc;
addpath('Functions\')

%folder with raw data files
datafolder = 'Data\';


%% import raw data

%if the data is in CDF files
T = importCDF(datafolder);

%make sampleIDs shorter and sort
T.sampleID = cellfun(@(x) upper(x(6:10)),T.sampleID,'UniformOutput',false);
T = sortrows(T,'sampleID');

%plot the raw data
figure; 
plot(T.RetentionTimeRaw(1,:),T.SignalRaw);
xlabel('Retention time (min)'); title('Raw chromatograms');
legend(T.sampleID)



%% Crop, baseline correction, normalization

retention_limits = [5 20];
T=cropAndInterpolate(T,retention_limits);

%%if you want to change the baseline parameters:
baseline_parameters.smoothness = 1e9;
baseline_parameters.asymmetry = 1e-4;
T = estimateBaseline(T,baseline_parameters);

figure; plot(T.RetentionTime(1,:),T.Baseline)

T=correctAndNormalize(T);

%plot the normalized chromatograms
figure; 
plot(T.RetentionTime(1,:),T.SignalNormalized);
xlabel('Retention time (min)'); title('Baseline corrected and normalized chromatograms');
legend(T.sampleID)


%% fit calibration model

% Import calibration data
caldata = readtable(strcat(datafolder,'calibration.xlsx'));

%fit calibration model 
calmodel = fitCalibrationModel(caldata.retention_time,log10(caldata.molecular_weight));


%plot normalized chromatograms together with calibration curve
figure;
xx = linspace(calmodel.MinTime-1,calmodel.MaxTime+3);
[~,yy,slope] = retentiontimeToMolarmass(xx,calmodel);
yyaxis left
plot(T.RetentionTime(1,:),T.SignalNormalized,'-');
ylabel('Normalized detector signal')
yyaxis right
semilogy(calmodel.xdata,10.^calmodel.ydata,'o','MarkerFaceColor','r')
hold on
semilogy(xx,10.^yy,'-r','linewidth',2)
ylabel('Molar mass')
xlabel('Retention time (min)'); title('Baseline corrected and normalized chromatograms');

% Add calibration model to table
T.CalModel = repmat(calmodel,height(T),1);


%% Calculate molar mass distributins and moments (Mw)

T = calculateMolweightdistr(T);

%plot
figure; semilogx(T.MolarMass(1,:),T.xM)
xlabel('Molar mass (g/mol)'); title('Differential log molecular weight distribution');
legend(T.sampleID)



%% Multiarea analysis

%convert retention times to molarmasses, for multiarea calculation
TimeLimits = [7.8 8.6 9.2 10 11]; %minutes
MolarMassLimits = retentiontimeToMolarmass(TimeLimits,calmodel);

%or define molar mass limits directly
MolarMassLimits = [100 300 600]


T = multiAreaAnalysis(T,MolarMassLimits);


figure; 
subplot(1,2,1); bar(T.Mw); set(gca,"XTickLabel",T.sampleID); ylabel('g/mol'); title('Average molecular weight')
subplot(1,2,2); bar(table2array(T(:,contains(T.Properties.VariableNames,'Fraction'))),'stacked'); ylabel('%')
legend(T.Properties.VariableDescriptions(contains(T.Properties.VariableNames,'Fraction')),'location','eastoutside')
set(gca,"XTickLabel",T.sampleID);





