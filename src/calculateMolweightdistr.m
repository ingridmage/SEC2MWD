function T=calculateMolweightdistr(T)
% T=calculateMolweightdistr(T)
% 
% INPUT
% T is a table that has to contain following variables:
%   RetentionTime           Retention times
%   SignalNormalized        Baseline corrected and area normalized detector signal
%   CalModel                struct generated by fitCalModel.m, with fields
%       -Type: (name of fit function)
%       -MinTime: minimum retention time for calibration data
%       -MaxTime: maximum retention time for calibration data
%       -p: fitted model coefficients
%       -extrapol_p1: coefficients for linear extrapolation below MinTime
%       -extrapol_p2: coefficients for linear extrapolation above MinTime
%              
%
% OUTPUT
% The same table with additional variables:
%
% M                 Molar mass
% logM              log10 of Molar Mass
% xM                differential molecular weight distribution
% wM                weight molecular weight distribution
% nM                number molecular weight distribution
% Mw                average weight molecular weight
% Mn                average number molecular weight
% PDI               Polydispersive index, a measure of heterogenity
% 
% Calculations are done as in references: 
% Gavrilov, M., & Monteiro, M. J. (2015). Derivation of the molecular weight distributions from size exclusion chromatography. European Polymer Journal, 65, 191–196. https://doi.org/10.1016/J.EURPOLYMJ.2014.11.018
% Shortt, D. W. (1993). Differential molecular weight distributions in high performance size exclusion chromatography. Journal of Liquid Chromatography, 16(16), 3371–3391. https://doi.org/10.1080/10826079308019695

if sum(contains(T.Properties.VariableNames,'RetentionTime'))==0; error('Input table must contain RetentionTime'); end
if sum(contains(T.Properties.VariableNames,'SignalNormalized'))==0; error('Input table must contain SignalNormalized'); end
if sum(contains(T.Properties.VariableNames,'CalModel'))==0; error('Input table must contain CalModel'); end


%remove potential "old" weigh distribution variables from T
T(:,strcmp(T.Properties.VariableNames,'logMolarMass'))=[];
T(:,strcmp(T.Properties.VariableNames,'MolarMass'))=[];
T(:,strcmp(T.Properties.VariableNames,'xM'))=[];
T(:,strcmp(T.Properties.VariableNames,'wM'))=[];
T(:,strcmp(T.Properties.VariableNames,'nM'))=[];
T(:,strcmp(T.Properties.VariableNames,'Mw'))=[];
T(:,strcmp(T.Properties.VariableNames,'Mn'))=[];
T(:,strcmp(T.Properties.VariableNames,'PDI'))=[];
T(:,strcmp(T.Properties.VariableNames,'FV'))=[];
T(:,strcmp(T.Properties.VariableNames,'WM'))=[];
T(:,strcmp(T.Properties.VariableNames,'slope'))=[];



for i=1:height(T)
    
   Res(i) = mwdistr(T.RetentionTime(i,:),T.SignalNormalized(i,:),T.CalModel(i));
   
    
end

Res = struct2table(Res);
T = [T Res];

%set units
T.Properties.VariableUnits{strcmp(T.Properties.VariableNames,'Mw')} = 'g/mol';
T.Properties.VariableUnits{strcmp(T.Properties.VariableNames,'Mn')} = 'g/mol';

%set descriptions
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'xM')} = 'differential log weight distribution';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'wM')} = 'differential weight distribution';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'PDI')} = 'Polydispersity index';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'FV')} = 'Cumulative weight fraction based on molecular weight';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'WM')} = 'Cumulative weight fraction based on retention time';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'Mw')} = 'Weight average molecular weight';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'Mn')} = 'Number average molecular weight';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'slope')} = 'Slope of column calibration curve';


end




function Results = mwdistr(Time,SignalNorm,CalModel) 

%transpose to column vectors if they are given as rows
if size(Time,1)==1; Time=Time'; end
if size(SignalNorm,1)==1; SignalNorm=SignalNorm'; end

%Check that V and H have same length
if length(Time)~=length(SignalNorm); error('Time and Signal must have same length'); end


%% Calculate M and slope from calibration curve

[M,logM,slope] = retentiontimeToMolarmass(Time,CalModel);

%% Calculate the cumulative distributions F(V) and W(M)
FV = cumtrapz(Time,SignalNorm);
WM = 1-FV;

%% calculate x(M), w(M) and n(M)
xM = -1./slope.*SignalNorm; 

% Flip logM, M and xM so that they are sorted accoring to ascending
%molecular weight. This is needed to obtain positive areas with the trapz
%function.
logM=flipud(logM);
M=flipud(M);
xM=flipud(xM);
WM = flipud(WM);

%calculate w(M) 
wM = xM./M*log10(exp(1));
Mw = trapz(M,wM.*M); %Average molecular weigh based on weight

filter = logM<0;
logM(filter)=nan;
M(filter)=nan;
xM(filter)=nan;
wM(filter)=nan;


%calculate n(M) 
nM = wM./M;
AnM = trapz(M(~isnan(wM)),nM(~isnan(wM)));

 

%calculate average molecular weights
Mn = 1/AnM; %Average moleculat weight based on number
PDI = Mw/Mn; %Polydispersive index, a measure of size heterogenity


Results.slope = slope';
Results.MolarMass = M';
Results.FV = FV';
Results.WM = WM';
Results.xM = xM';
Results.wM = wM'; 
Results.nM = nM'; 
Results.Mw = Mw;
Results.Mn = Mn;
Results.PDI = PDI;

end