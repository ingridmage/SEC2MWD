function T = correctAndNormalize(T)
% T=correctAndNormalize(T)
%
% Subtracts baseline and normalizes to unit area under curve
%
% INPUT
% T                         Table that has to contain following variables (Table can be obtained my importCDF);
%                           - SignalCropped                  Detector signal, cropped to area of interest
% OUTPUT
% The same table with additional variables:
%
% SignalNormalized          Baseline corrected and normalized detector signal (area = 1)
% TotalArea                 Area under the curve after cropping, used for normalization


if sum(contains(T.Properties.VariableNames,'SignalCropped'))==0; error('Input table must contain SignalCropped'); end
if sum(contains(T.Properties.VariableNames,'Baseline'))==0; error('Input table must contain Baseline'); end

signal = T.SignalCropped-T.Baseline;
signal(signal<0)=0; %set negative values to zero

T.TotalArea = trapz(T.RetentionTime(1,:)',signal')';
T.SignalNormalized = signal./repmat(T.TotalArea,1,size(signal,2));



%variable desriptions
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'SignalNormalized')} = 'Baseline corrected and normalized detector signal';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'TotalArea')} = 'Area under baseline corrected chromatogram, within retention time limits';



