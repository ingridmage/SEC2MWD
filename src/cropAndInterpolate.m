function T=cropAndInterpolate(T,retention_limits,resolution)
% T=cropAndInterpolate(T,retention_limits,resolution)
%
% This function crops the chromatogram and changes the resolution
%
% INPUT
% T                         Table that has to contain following variables (Table can be obtained my importCDF);
%                           - RetentionTimeRaw           Retention times, one vector per table row
%                           - SignalRaw                  Detector signal, same dimension as RetentionTimeRaw
%
% retention_limits           optional 2-element vector. Default is [5 20], i.e. crop to 5-20 minutes
%
% resolution                optional. Deault is 1800
%
% OUTPUT
% The same table with additional variables:
%
% RetentionTime             cropped to 5-20 minutes
% SignalCropped             Cropped detector signal, interpolated to the given resolution

if ~exist('retention_limits', 'var') | isempty(retention_limits)
        retention_limits = [5 20];
end


if ~exist('resolution', 'var') | isempty(resolution)
        resolution = 1800;
end

%check if table contains the necessary variables
if sum(contains(T.Properties.VariableNames,'RetentionTimeRaw'))==0; error('Input table must contain RetentionTimeRaw'); end
if sum(contains(T.Properties.VariableNames,'SignalRaw'))==0; error('Input table must contain SignalRaw'); end


%remove potential "old" result variables from T
T(:,strcmp(T.Properties.VariableNames,'RetentionTime'))=[];
T(:,strcmp(T.Properties.VariableNames,'SignalCropped'))=[];

%initialize new variables
T.RetentionTime = nan(height(T),resolution);
T.SignalCropped = nan(height(T),resolution);


for i=1:height(T)
    
    V=T.RetentionTimeRaw(i,:);
    H=T.SignalRaw(i,:);
    
    %filter for cropping
    filter = (V>retention_limits(1) & V<retention_limits(2));

    %interpolate to resolution
    RetTime = linspace(retention_limits(1),retention_limits(2),resolution);
    Signal = interp1(V(filter),H(filter),RetTime,'linear','extrap');
    
    
    T.RetentionTime(i,:)= RetTime;
    T.SignalCropped(i,:)=Signal;
    
end


%set units and descriptions
T.Properties.VariableUnits{strcmp(T.Properties.VariableNames,'RetentionTimeRaw')} = 'minutes';
T.Properties.VariableUnits{strcmp(T.Properties.VariableNames,'RetentionTime')} = 'minutes';

T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'RetentionTime')} = 'Retention time';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'SignalCropped')} = 'Cropped and interpolated raw detector signal';

end


