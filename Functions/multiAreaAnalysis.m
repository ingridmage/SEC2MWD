function T = multiAreaAnalysis(T,molarmass_limits,plots)
% T = multiarea(T,molarmass_limits,plots)
%
% Limits must be in molar masses (not retention time) since only molar
% masses are comparable if samples are based on different column
% calibrations.
%
%
% INPUT
%
% T                     Table 
% molarmass_limits      Limits between fractions (g/mol). Do not include end points
% plots                 =1 outputs plot of calibration curve and derivative
%
% OUTPUT
% T                     Table with new "fraction" columns

if nargin == 1
    error('Need to input both Table and molarmass_limits')
elseif nargin == 2
    plots = 1; %default plot
end

molarmass_limits = sort(molarmass_limits);

%remove "old" area variables
T(:,contains(T.Properties.VariableNames,'Fraction'))=[];

%find indices for fraction limits
idx =[];
for i=1:length(molarmass_limits)
    [~,idx(i)] = min(abs(T.logMolarMass(1,:)-log10(molarmass_limits(i))));
end

idx = ([1 idx length(T.WM(i,:))]); %include min and max


%lookup areas from cumulative weight distribution
for i = 1:height(T)    
    A(i,:) = diff([T.WM(i,idx)*100]);
end


%add to table
for i=1:size(A,2)
   varname = strcat('Fraction',num2str(size(A,2)-i+1));
   T.(varname) = A(:,i);
   T.Properties.VariableUnits{strcmp(T.Properties.VariableNames,varname)} = '%';
   
   if i==1 %lowest fraction
       ss = ['<' num2str(round(molarmass_limits(i))) ' g/mol'];
   T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,varname)} = ss;
   elseif i==size(A,2) %highest fraction
          ss = ['>' num2str(round(molarmass_limits(end))) ' g/mol'];
   T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,varname)} = ss;
   else 
       ss = [num2str(round(molarmass_limits(i-1))) '-' num2str(round(molarmass_limits(i))) ' g/mol'];
          T.Properties.VariableDescriptions{strcmp(T.Properties.VariableNames,varname)} = ss;
   end
       
      
    
end


if plots ==1
    
    figure;
    semilogx(T.MolarMass(1,:),T.xM,'k')
    xline(molarmass_limits,'r')
    xlabel('log_{10} Molar Mass')
    ylabel('differential log weight distribution')
    
end






