function calmodel = fitCalibrationModel(x,y,fitfunction,plotresults)
% calmodel = fitCalibrationModel(x,y,fitfunction,plotresults)
%
% fits calibration model, with linear extrapolation outside the calibration
% data range
%
% INPUT
% x             vector of retention times
% y             vector of log10 molecular weights (same size as x)
% fitfunction   fitfunction(s), scalar or vector:
%                   0 = all ,
%                   1 = linear,
%                   2 = Poly3
%                   3 = PSS Poly3 (default)
% plotresults =1 (default) creates a plot of the fitted curve
%
%
% OUTPUT
% calmodel    struct with fields:
%              -p:          model parameters. It is a vector with length (polynomial order +1).
%                           First element is highst polynomial order, last is intercept.
%                           Example: CalParam = [a b c d] means log(M) = a*V^3 + b*V^2 + c*V + d
%              -MinTime:     Minimum Retention time for calibration samples
%              -MaxTime:     Maximum Retention time for calibration samples
%              -extrapol_p1  Coefficients for linear interpolation below minimum calibration sample. First element is slope, second is intercept
%              -extrapol_p2  Coefficients for linear interpolation above maximum calibration sample. First element is slope, second is intercept

%check inputs and set defaults
if ~exist('x') | ~exist('y'); error('Both x and y must be provided'); end

if size(x,1) == size(y,2); x=x'; end

if size(x)~=size(y); error('x and y must have same size'); end

if ~exist('fitfunction','var'); fitfunction = 3; end
if ~exist('plotresults','var'); plotresults = 1; end

%remove samples with nan
id = isnan(x) | isnan(y);
if sum(id)>0
    x(id)=[]; y(id)=[];
    disp(['Removing ' num2str(sum(id)) ' calibration sample(s)  due to missing values' ])
end




%fit calibration model(s)
minTime=min(x);
maxTime=max(x);

calmodel = struct();
i=1;

if sum(fitfunction == 0)==1 | sum(fitfunction ==1)>0 %fit linear model

    calmodel(i).Type = 'Linear';
    calmodel(i).MinTime = minTime;
    calmodel(i).MaxTime = maxTime;

    calmodel(i).p = polyfit(x,y,1);

    calmodel(i).xdata = x;
    calmodel(i).ydata = y;

    i=i+1;

end

if sum(fitfunction == 0)==1 | sum(fitfunction ==2)>0 %fit "Poly3", regular 3rd order polynomial

    calmodel(i).Type = 'Poly3';
    calmodel(i).MinTime = minTime;
    calmodel(i).MaxTime = maxTime;

    calmodel(i).p = polyfit(x,y,3);

    calmodel(i).xdata = x;
    calmodel(i).ydata = y;

    i=i+1;
end

if sum(fitfunction == 0)==1 | sum(fitfunction ==3)>0 %fit "PSS Poly3", average between linear and 3rd order polynomial

    calmodel(i).Type = 'PSS Poly 3';
    calmodel(i).MinTime = minTime;
    calmodel(i).MaxTime = maxTime;

    p1=polyfit(x,y,1);
    p3=polyfit(x,y,3);



    calmodel(i).p = mean([0 0 p1;p3]);

    calmodel(i).xdata = x;
    calmodel(i).ydata = y;
    i=i+1;
end


%make linear extrapolations
%define dummy Time-vector for calculation of extrapolation
dummyTime = [minTime-0.01 minTime minTime+0.01 maxTime-0.01 maxTime maxTime+0.01];
for i=1:length(calmodel)

    %fit dummyTime values
    dummylogM = polyval(calmodel(i).p,dummyTime);


    %make linear interpolation below minT, using 3 points
    p1=polyfit(dummyTime(1:3),dummylogM(1:3),1);
    calmodel(i).extrapol_p1 = p1;

    %make linear interpolation above maxT, using 3 points
    p2=polyfit(dummyTime(4:6),dummylogM(4:6),1);
    calmodel(i).extrapol_p2 = p2;


end



if plotresults ==1
    figure;

    xx = linspace(calmodel(1).MinTime-1,calmodel(1).MaxTime+5);
    for i=1:length(calmodel)
        yy=polyval(calmodel(i).p,xx);
        yyder = polyval(polyder(calmodel(i).p),xx);


        %fit extrapolations
        yy(xx<minTime) = polyval(calmodel(i).extrapol_p1,xx(xx<minTime));
        yyder(xx<minTime) = polyval(polyder(calmodel(i).extrapol_p1),xx(xx<minTime));

        yy(xx>maxTime) = polyval(calmodel(i).extrapol_p2,xx(xx>maxTime));
        yyder(xx>maxTime) = polyval(polyder(calmodel(i).extrapol_p2),xx(xx>maxTime));

        %plot
        subplot(1,2,1)
        plot(xx,yy)
        hold on
        subplot(1,2,2)
        plot(xx,yyder)
        hold on
    end
    subplot(1,2,1)
    legend(calmodel.Type,'AutoUpdate','off')
    scatter(x,y,'k','filled');

    xlabel('Retention time')
    ylabel('log10 molar mass')

    subplot(1,2,2)
    xlabel('Retention time')
    ylabel('derivative')
end


