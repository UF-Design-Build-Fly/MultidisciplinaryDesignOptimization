function SaveTopPlanesCSV(topPlanes, filename)
% SaveTopPlanesCSV  –  writes a summary CSV for quick review/sharing.

    if nargin < 2 || isempty(filename), filename = 'topPlanes_2026.csv'; end
    K = numel(topPlanes);
    if K==0, warning('No planes to export.'); return; end

    % Preallocate
    Span_m = zeros(K,1); AR = zeros(K,1); RAC = zeros(K,1);
    Motor = strings(K,1); Cells = zeros(K,1); PropD_in = zeros(K,1); PropP_in = zeros(K,1);
    Wh = zeros(K,1);
    P = zeros(K,1); C = zeros(K,1); Banner_in = zeros(K,1); BannerH_in = zeros(K,1);
    V2 = zeros(K,1); V3 = zeros(K,1); Laps2 = zeros(K,1); Laps3 = zeros(K,1);
    M2 = zeros(K,1); M3 = zeros(K,1); GM = zeros(K,1); Score = zeros(K,1);
    mEmpty = zeros(K,1); m2Gross = zeros(K,1); m3Gross = zeros(K,1);

    for i = 1:K
        p = topPlanes(i);
        Span_m(i)  = safe(p, 'wing', 'span');
        AR(i)      = safe(p, 'wing', 'aspectRatio');

        WS_in = Span_m(i)/0.0254; WS_in_meas = round(WS_in); WS_ft_rec = WS_in_meas/12;
        RAC(i) = max(0.90, 0.05*WS_ft_rec + 0.75);

        Motor(i)   = trystr(p, 'powerSystem', 'motorName');
        Cells(i)   = safe(p, 'powerSystem', 'cells');
        PropD_in(i)= safe(p, 'powerSystem', 'propDiameter');
        PropP_in(i)= safe(p, 'powerSystem', 'propPitch');
        Wh(i)      = safe(p, 'powerSystem', 'batteryCapacity');

        P(i)       = safe(p, 'performance', 'passengers');
        C(i)       = safe(p, 'performance', 'cargo');
        Banner_in(i)= safe(p, 'performance', 'bannerLength_in');
        BannerH_in(i)= safe(p, 'performance', 'bannerHeight_in');

        V2(i)      = safe(p, 'performance', 'velocity2');
        V3(i)      = safe(p, 'performance', 'velocity3');
        Laps2(i)   = safe(p, 'performance', 'laps2');
        Laps3(i)   = safe(p, 'performance', 'laps3');

        M2(i)      = safe(p, 'performance', 'M2');
        M3(i)      = safe(p, 'performance', 'M3');
        GM(i)      = safe(p, 'performance', 'GM');
        Score(i)   = safe(p, 'performance', 'CompetitionScore');

        mEmpty(i)  = safe(p, 'performance', 'totalEmptyWeight');
        m2Gross(i) = safe(p, 'performance', 'totalWeight2');
        m3Gross(i) = safe(p, 'performance', 'totalWeight3');
    end

    T = table((1:K)', Span_m, AR, RAC, Motor, Cells, PropD_in, PropP_in, Wh, ...
              P, C, Banner_in, BannerH_in, V2, V3, Laps2, Laps3, M2, M3, GM, ...
              mEmpty, m2Gross, m3Gross, Score, ...
              'VariableNames', {'Rank','Span_m','AR','RAC','Motor','Cells','PropD_in','PropP_in','Wh', ...
                                'P','C','Banner_in','BannerH_in','V2','V3','Laps2','Laps3','M2','M3','GM', ...
                                'EmptyMass_kg','M2_Mass_kg','M3_Mass_kg','CompetitionScore'});

    writetable(T, filename);
    fprintf('CSV written: %s\n', filename);

    function val = safe(s, f1, f2)
        try
            val = s.(f1).(f2);
            if ~isnumeric(val), val = double(val); end
        catch
            val = NaN;
        end
    end
    function str = trystr(s, f1, f2)
        try
            str = string(s.(f1).(f2));
        catch
            str = "";
        end
    end
end
