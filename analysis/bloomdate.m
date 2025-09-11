function A = bloomdate(In, var, opt)

arguments
    In {mustBeA(In, "timetable")}
    var {mustBeTextScalar}
    opt.thresh (1,1) {mustBeNumeric} =1.05
    opt.doylimit (1,1) {mustBeNumeric, mustBeLessThan(opt.doylimit,367), mustBeNonnegative} =250
    opt.plotflag (1,1) {mustBeNumericOrLogical} =false
    opt.fracmax (1,1) {mustBeNumeric, mustBeLessThan(opt.fracmax,1.0), mustBeNonnegative} =0.5
end

Tmp = timetable2table(In);

% Separate "years" by low points

yrfrac = 0.8;
[pklo, tlo] = findpeaks(-Tmp.(var), Tmp.Time, 'minpeakdistance', years(yrfrac), 'minpeakwidth', days(20), 'annotate', 'extents');

% if opt.plotflag
%     figure;
%     plot(In, var);
%     hold on;
%     plot(tlo, -pklo, '^');
% end

% yrlim = year([min(In.Time), max(In.Time)]);
% mask = tlo > datetime(yrlim(1),7,1) & tlo < datetime(yrlim(2),7,1);
% tlo = tlo(mask); % throwing out end years, partial in some models
% pklo = pklo(mask);

% Add table variables indicating "year" bin and day of "year"

Tmp.yrbin = discretize(Tmp.Time, tlo);
Tmp = Tmp(~isnan(Tmp.yrbin),:);

yrnom = year(splitapply(@mode, dateshift(Tmp.Time, 'start', 'year'), Tmp.yrbin));

Tmp.doy = days(Tmp.Time - datetime(yrnom(Tmp.yrbin),1,1));
Tmp.yrnom = yrnom(Tmp.yrbin);

A = rowfun(@(x,y) bloomstats(x,y,opt.thresh,opt.doylimit,opt.fracmax), Tmp, "GroupingVariable","yrbin", ...
           "InputVariables",["doy" var], ...
           "OutputVariableNames",["max_val","max_doy", ...
                                  "bloom_val","bloom_doy",...
                                  "min_val","min_doy", ...
                                  "half_val", "half_doy", ...
                                  "maxx_val", "maxx_doy"]);

B = rowfun(@(x) mean(x,1), Tmp, "GroupingVariable","yrbin", "InputVariables", "yrnom", "OutputVariableName", "yrnom");

A.yrnom = B.yrnom;

if opt.plotflag
    n = max(Tmp.yrbin);
    ncol = 4;
    nrow = ceil(n/ncol);

    h = plotgrid('size', [nrow ncol], 'sp', 0.01);
    for ii = 1:max(Tmp.yrbin)
        Plt = Tmp(Tmp.yrbin==ii,:);

        plot(h.ax(ii), Plt.doy, Plt.(var), '-', ...
            A.max_doy(ii), A.max_val(ii), 'o', ...
            A.bloom_doy(ii), A.bloom_val(ii), '*', ...
            A.half_doy(ii), A.half_val(ii), 'x');
    end
    set(h.fig, 'color', 'w');
    set(h.ax, 'box', 'off', 'ycolor', [1 1 1]*0.8, 'xcolor', [1 1 1]*0.8, 'xlim', [-50 400], 'clipping', 'off');
    set(h.ax(:,2:end), 'yticklabel', ''); set(h.ax(1:end-1,:), 'xticklabel', '');
    labelaxes(h.ax(1:n), compose("%d", A.yrnom), 'northwest');
end

end

% Subfunction: calculate peak and "start" date for bloom

function [ymax, tmax, ybloom, tbloom, ymin, tmin, yhalf, thalf, ymaxx, tmaxx] = bloomstats(times, vals, thresh, doylimit, fracmax)
%BLOOMSTATS Identifies key points in spring bloom
%
% Identifies time and value associated with a) peak (spring) value, b)
% bloom start, i.e. time when values rise above a threshold relative to
% median, and c) start, i.e. time when rate of change becomes positive

    % Max

    valstmp = vals;
    valstmp(times>doylimit & ~isnan(vals)) = -Inf;
    [ymax,imax] = max(valstmp);

    if ~isnan(ymax)
        tmax = times(imax);
    else
        tmax = NaN;
    end

    [ymaxx,imaxx] = max(vals);

    if ~isnan(ymaxx)
        tmaxx = times(imaxx);
    else
        tmaxx = NaN;
    end

    
    % Bloom (crosses threshold)

    threshval = median(vals,'omitnan').*thresh;

    isunder = vals < threshval & times < times(imax);
    if any(isunder)
        ibloom = find(isunder, 1, 'last');
        tbloom = times(ibloom+1);
        ybloom = vals(ibloom+1);
    else
        tbloom = NaN;
        ybloom = NaN;
    end

    % Min
    
    if all(isnan(vals))
        ymin = NaN;
        tmin = NaN;
    else
        ymin = vals(1);
        tmin = times(1);
    end

    % Half (reaches half peak value)

    threshval = ymin + (ymax-ymin)*fracmax;

    isunder = vals < threshval & times < times(imax);
    if any(isunder)
        ihalf = find(isunder, 1, 'last');
        thalf = times(ihalf+1);
        yhalf = vals(ihalf+1);
    else
        thalf = NaN;
        yhalf = NaN;
    end


end


