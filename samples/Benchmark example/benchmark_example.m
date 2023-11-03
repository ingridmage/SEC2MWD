%% this example is taken from
% 
% Gavrilov, M., & Monteiro, M. J. (2015). 
% Derivation of the molecular weight distributions from size exclusion chromatography. 
% European Polymer Journal, 65, 191–196. 
% https://doi.org/10.1016/J.EURPOLYMJ.2014.11.018
% 
% where the data is given as supporting material 

%% 

clear; close all; clc
addpath('..\..\src')

%% import data

Ve = readmatrix("benchmark_example.xlsx",'Range','I2:I326');
hv = readmatrix("benchmark_example.xlsx",'Range','L2:L326');

T=cell2table({'Sample' Ve' hv'},'VariableNames', {'sampleID','RetentionTimeRaw','SignalRaw'})

%plot the raw data
figure; 
plot(T.RetentionTimeRaw(1,:),T.SignalRaw);
xlabel('Retention time (min)'); title('Raw chromatograms');
legend(T.sampleID)



%% Crop, baseline correction, normalization

%the full range is used
retention_limits = [min(T.RetentionTimeRaw) max(T.RetentionTimeRaw)];
T=cropAndInterpolate(T,retention_limits);

% No baseline correction was done, so the estiamted baseline is set to zero 
T = estimateBaseline(T);
T.Baseline = T.Baseline*0;

T=correctAndNormalize(T);

%plot the normalized chromatograms
figure; 
plot(T.RetentionTime(1,:),T.SignalNormalized);
xlabel('Retention time (min)'); title('Baseline corrected and normalized chromatograms');
legend(T.sampleID)


%% fit calibration model


calmodel.Type = 'Poly5';
calmodel.MinTime = min(T.RetentionTime(1,:));
calmodel.MaxTime = max(T.RetentionTime(1,:));
calmodel.p = readmatrix("benchmark_example.xlsx",'Range','D2:D7')';
calmodel.extrapol_p1 = [];
calmodel.extrapol_p2 = [];





%plot normalized chromatograms together with calibration curve
figure;
xx = linspace(calmodel.MinTime,calmodel.MaxTime);
[~,yy,slope] = retentiontimeToMolarmass(xx,calmodel);
yyaxis left
plot(T.RetentionTime(1,:),T.SignalNormalized,'-');
ylabel('Normalized detector signal')
yyaxis right
semilogy(xx,10.^yy,'-r','linewidth',2)
ylabel('Molar mass')
xlabel('Retention time (min)'); title('Baseline corrected and normalized chromatograms');

% Add calibration model to table
T.CalModel = repmat(calmodel,height(T),1);


%% Calculate molar mass distributins and moments (Mw)

T = calculateMolweightdistr(T);

%plot
f=figure; 


subplot(3,2,1)
xx = linspace(calmodel.MinTime,calmodel.MaxTime);
[~,yy,slope] = retentiontimeToMolarmass(xx,calmodel);
yyaxis left
plot(T.RetentionTime(1,:),T.SignalNormalized,'LineWidth',2);
ylabel('chromatogram')
yyaxis right
semilogy(xx,10.^yy,'-r','linewidth',2)
ylabel('Molar mass')
xlabel('Elution volume'); 

subplot(3,2,2)
semilogx(T.MolarMass(1,:),T.xM,'LineWidth',2)
xlabel('Molar mass (g/mol)'); title('xM');

%plot
subplot(3,2,3)
semilogx(T.MolarMass(1,:),T.wM,'LineWidth',2)
xlabel('Molar mass (g/mol)'); title('wM');

%plot
subplot(3,2,4)
 semilogx(T.MolarMass(1,:),T.nM,'LineWidth',2)
xlabel('Molar mass (g/mol)'); title('nM');

tmp = T(:,{'Mw' 'Mn' 'PDI'});
tmp.Mn = round(tmp.Mn);
tmp.PDI = round(tmp.PDI,2)

TString = evalc('disp(tmp)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');

annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','Position',[0.3500    0.0500    0.3339    0.2000]);









