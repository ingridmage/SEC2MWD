function T = importCDF(folder)
% Imports cdf files from Thermo Fisher SEC
% 
% INPUT: folder (string)  full path of folder containing the cdf files
% OUTPUT T      (table)   Table with sample info and raw data for all cdf-files in folder
%
%
% file format is NetCDF. Type 'help netcdf' for documentation on file type.


%make sure the path name ends with  \
if ~strcmp(folder(end),'\')
    folder = strcat(folder,'\');
end
    

%make list of all cdf files in folder
fileList = dir(strcat(folder,'*.cdf'));
nFiles = length(fileList);


%initialize struct
T=struct();


if nFiles>0
    disp(['Importing raw data from ' num2str(nFiles) ' cdf-files to MATLAB table'])

%read sample attributes
for i=1:nFiles
    
        if round(i/10)==i/10; disp(['Importing file ' num2str(i) ' of ' num2str(nFiles)]); end

    
    file = fullfile(fileList(i).folder,fileList(i).name);
    
    
    %read file info
    info = ncinfo(file);
    
   %sample name may be found in different fields
    if any(strcmpi('experiment_title', {info.Attributes.Name}))
        T(i).sampleID = ncreadatt(file, '/', 'experiment_title');
    end

    if any(strcmpi('sample_name', {info.Attributes.Name}))
        T(i).sampleName = ncreadatt(file, '/', 'sample_name');
    end

    
    %file name
    T(i).FileName = strrep(fileList(i).name,'.cdf','');

    
    %read retention times and detector response
    if isfield(info, 'Variables')
        
        % Check for time values
        if any(strcmpi('raw_data_retention', {info.Variables.Name}))
            T(i).RetentionTimeRaw=ncread(file, 'raw_data_retention')' ./60;
        else %retention times are not given explicitely

            start = ncread(file,'actual_delay_time');
            stop = ncread(file,'actual_run_time_length');
            int = ncread(file,'actual_sampling_interval');
            T(i).RetentionTimeRaw=(start:int:stop)./60;


        end
        
        % Check for total intensity values
        if any(strcmpi('ordinate_values', {info.Variables.Name}))
            T(i).SignalRaw= ncread(file, 'ordinate_values')';
        end
        
        if length( T(i).RetentionTimeRaw)~=length(T(i).SignalRaw)
            error(['Length of retention time cevort not equal to length of signal vector for file ' file])
        end
        
    end
    
    if any(strcmpi('injection_date_time_stamp', {info.Attributes.Name}))
        tmp = ncreadatt(file, '/', 'injection_date_time_stamp');
        if contains(tmp,'+'); tmp = strsplit(tmp,'+'); tmp = tmp{1}; end
        T(i).InjectionTime = datetime(tmp,'InputFormat','yyyyMMddHHmmss');
    end
    
    %some extra info
    
    if any(strcmpi('dataset_origin', {info.Attributes.Name}))
        T(i).DataOrigin = ncreadatt(file, '/', 'dataset_origin');
    end
    
    
    if any(strcmpi('detector_name', {info.Attributes.Name}))
        T(i).DetectorName = ncreadatt(file, '/', 'detector_name');
    end
    
    if any(strcmpi('sample_injection_volume', {info.Attributes.Name}))
        T(i).SampleInjectionVolume = ncreadatt(file, '/', 'sample_injection_volume');
    end
    
    
    
end

%interpolate to same retention times.One second resolution
maxtime = max(cell2mat(arrayfun(@(x)(max(x.RetentionTimeRaw)), T, 'UniformOutput', false)));
RetTime = 0:1/60:maxtime;

for i=1:length(T)
    
     T(i).SignalRaw  = interp1(T(i).RetentionTimeRaw,T(i).SignalRaw,RetTime);
    T(i).RetentionTimeRaw = RetTime; 
end



T=struct2table(T);

T.Properties.VariableUnits{contains(T.Properties.VariableNames,'RetentionTime')} = 'minutes';
T.Properties.VariableDescriptions{contains(T.Properties.VariableNames,'RetentionTime')} = 'Retention Time';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'SignalRaw')} = 'Raw detector signal';
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'InjectionTime')} = 'date and time of injection';
T.Properties.VariableUnits{strcmp(T.Properties.VariableNames,'SampleInjectionVolume')} = 'ml';

end

if nFiles == 0; disp('Found no CDF files in folder'); end

end

