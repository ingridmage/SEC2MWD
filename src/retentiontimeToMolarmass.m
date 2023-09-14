function [MolarMass,LogMolarMass,slope] = retentiontimeToMolarmass(Time,calmodel)
% [MolarMass,LogMolarMass,slope] = retentiontimeToMolarmass(Time,calmodel)
%
% INPUT
% Time          Vector of retention times
% calmodel      Struct created by function 'fitcalmodel'
%
% OUTPUT
% MolarMass     MolarMass corresponding to Time
% LogMolarMass  Log10 MolarMass
% slope         The derivative of the calibration curve at every Time




if strcmp(calmodel.Type,'Linear') | strcmp(calmodel.Type,'Poly3')| strcmp(calmodel.Type,'PSS Poly 3')
    LogMolarMass=polyval(calmodel.p,Time);
            slope = polyval(polyder(calmodel.p),Time);



elseif strcmp(calmodel.Type,'Exponential')
    LogMolarMass=calmodel.p(1) + calmodel.p(2) * exp(-Time *calmodel.p(3));
            Time2=[Time Time(end)+min(diff(Time))];
            slope=diff(calmodel.p(1) + calmodel.p(2) * exp(-xx2 *calmodel.p(3)))./diff(xx2);

end


%fit extrapolations
LogMolarMass(Time<calmodel.MinTime) = polyval(calmodel.extrapol_p1,Time(Time<calmodel.MinTime));
LogMolarMass(Time>calmodel.MaxTime) = polyval(calmodel.extrapol_p2,Time(Time>calmodel.MaxTime));
        
slope(Time<calmodel.MinTime) = polyval(polyder(calmodel.extrapol_p1),Time(Time<calmodel.MinTime));
slope(Time>calmodel.MaxTime) = polyval(polyder(calmodel.extrapol_p2),Time(Time>calmodel.MaxTime));



MolarMass=10.^LogMolarMass;