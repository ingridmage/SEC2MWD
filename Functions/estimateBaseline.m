function T = estimateBaseline(T,parameters)
% T = estimateBaseline(T,parameters)
%
% Baseline estimation using asymmetric least squares
% Reference: P.H.C. Eilers, Analytical Chemistry, 75 (2003) 3631
%
% This code is a modified verion of baseline correction method in the
% Chromatograpy toolbox for MATLAB:
% https://github.com/chemplexity/chromatography
% https://se.mathworks.com/matlabcentral/fileexchange/47696-chromatography-toolbox
%
%
% INPUT
% T                         Table that has to contain following variables (Table can be obtained my importCDF);
%                           - RetentionTime           Retention times, one vector per table row
%                           - SignalCropped           Detector signal cropped to range of interest
%
% baseline_parameters       optional struct with parameters for baseline estimation. Fields:
%                               -smoothness (default = 1e7)
%                               -asymmetry (default 1e-9)
%
% OUTPUT
% The same table with additional variable:
%
% Baseline                  Estimated baseline

if ~exist('baseline_parameters','var')
    parameters.smoothness = 1e7;
    parameters.asymmetry = 1e-9;
end



if sum(contains(T.Properties.VariableNames,'RetentionTime'))==0; error('Input table must contain RetentionTime'); end
if sum(contains(T.Properties.VariableNames,'SignalCropped'))==0; error('Input table must contain SignalCropped'); end


%remove potential "old" result variables from T
T(:,strcmp(T.Properties.VariableNames,'Baseline'))=[];


%estimate and subtract baseline. Using absolute value of signal to
%prevent negative peaks (artifacts) to affect the baseline estimation
T.Baseline = nan(height(T),size(T.SignalCropped,2));
for i=1:height(T)
    bl = assymetric_leastsquares(abs(T.SignalCropped(i,:)'), parameters);
    T.Baseline(i,:) = bl';
end
T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,'Baseline')} = 'Estimated baseline';







function baseline = assymetric_leastsquares(y, parameters)
% INPUT
%   y              intensity values (vector)
%   parameters     struct with parameters:
%                               -smoothness (default: 1E7, suggested range: 1E3 to 1E9)
%                               -asymmetry (default: 1E-9, suggested range:  1E-3 to 1E-9)
%
%   bl                estimated baseline


% Check input
if nargin < 1
    error('Not enough input arguments.');

elseif ~isnumeric(y)
    error('Undefined input arguments of type ''y''.');
end

if nargin<2
elseif ~isfield(parameters,'smoothness')
    error('"parameters" should be a struct with fields "smoothness" and "asymmetry"');
elseif ~isfield(parameters,'asymmetry')
    error('"parameters" should be a struct with fields "smoothness" and "asymmetry"');
end

smoothness = parameters.smoothness;
asymmetry = parameters.asymmetry;

% Check data precision
if ~isa(y, 'double')
    y = double(y);
end

% Check data length
if all(size(y) <= 3)
    error('Insufficient number of points.');
end

% Check for negative values
if any(min(y) < 0)

    % Determine offset
    offset = min(y);
    offset(offset > 0) = 0;
    offset(offset < 0) = abs(offset(offset < 0));

else
    offset = zeros(1, length(y(1,:)));
end

% Variables
rows = length(y(:,1));
index = 1:rows;

% Pre-allocate memory
baseline = zeros(size(y));
weights = ones(rows, 1);

w = spdiags(weights, 0, rows, rows);

% Variables
d = diff(speye(rows), 2);
d = smoothness * (d' * d);

% Calculate baseline
for i = 1:length(y(1,:))

    % Check offset
    if offset(i) ~= 0
        y(:,i) = y(:,i) + offset(i);
    end

    % Check values
    if nnz(y(:,i)) == 0
        continue
    end

    % Pre-allocate memory
    b = zeros(rows,1);

    % Number of iterations
    for j = 1:10

        [~,pd] = chol(w + d);
        if pd~=0
           % warning('Problems with baseline estimation')
            continue
        end

        % Cholesky factorization
        w = chol(w + d);

        % Left matrix divide, multiply matrices
        b = w \ (w' \ (weights .* y(:,i)));

        % Determine weights
        weights = asymmetry * (y(:,i) > b) + (1 - asymmetry) * (y(:,i) < b);

        % Reset sparse matrix
        w = sparse(index, index, weights);
    end

    % Remove negative values
    b(b<0) = 0;

    % Remove offset
    if offset(i) ~= 0
        baseline(:,i) = b - offset(i);

    elseif offset(i) == 0
        baseline(:,i) = b;
    end

    % Reset variables
    weights = ones(rows, 1);
end


