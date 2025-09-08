function DrawPayloadPlan(ax, fitInfo, P, C)
% DrawPayloadPlan  —  schematic top-down layout of pucks + ducks.
% Uses the partitioned lengths from CanFitPayload's info struct.
    DUCK_W = 2.3; DUCK_L = 2.5;
    PUCK_D = 3.0;

    Lbay = fitInfo.bay.length_in;
    Wbay = fitInfo.bay.width_in;

    % Zones along length (pucks, then ducks)
    Lp = fitInfo.puck.len_req_in;
    Ld = fitInfo.duck.len_req_in;

    if ~isfinite(Lp), Lp = 0; end
    if ~isfinite(Ld), Ld = 0; end
    Ltot = min(Lbay, Lp + Ld);

    % Draw bay
    rectangle(ax, 'Position',[0 0 Lbay Wbay], 'EdgeColor',[0 0 0], 'LineWidth',1.5);
    hold(ax,'on');

    % Draw puck zone
    rectangle(ax, 'Position',[0 0 Lp Wbay], 'FaceColor',[0.90 0.95 1.00], 'EdgeColor',[0.2 0.3 0.7]);
    text(ax, Lp/2, Wbay+0.2, sprintf('Pucks (C=%d)', C), 'HorizontalAlignment','center');

    % Draw duck zone
    rectangle(ax, 'Position',[Lp 0 Ld Wbay], 'FaceColor',[1.00 0.95 0.90], 'EdgeColor',[0.7 0.4 0.2]);
    text(ax, Lp + Ld/2, Wbay+0.2, sprintf('Ducks (P=%d)', P), 'HorizontalAlignment','center');

    % Grid the pucks (approx)
    cap_puck_row = floor(Wbay / PUCK_D);
    n_rows = max(1, cap_puck_row);
    n_cols = ceil(max(0,C) / n_rows);
    px = 0; py = 0;
    for i=1:min(C, n_rows*n_cols)
        col = ceil(i/n_rows); row = i - (col-1)*n_rows;
        x = px + (col-1)*PUCK_D; y = py + (row-1)*PUCK_D;
        if x + PUCK_D <= Lp && y + PUCK_D <= Wbay
            rectangle(ax, 'Position',[x y PUCK_D PUCK_D], 'Curvature',1, ...
                      'FaceColor',[0.75 0.85 1.00], 'EdgeColor',[0.2 0.3 0.7]);
        end
    end

    % Grid the ducks (approx, upright)
    cap_duck_row = floor(Wbay / DUCK_W);
    n_rows = max(1, cap_duck_row);
    n_cols = ceil(max(0,P) / n_rows);
    dx0 = Lp;  % start at duck zone
    for i=1:min(P, n_rows*n_cols)
        col = ceil(i/n_rows); row = i - (col-1)*n_rows;
        x = dx0 + (col-1)*DUCK_L; y = (row-1)*DUCK_W;
        if x + DUCK_L <= Lp+Ld && y + DUCK_W <= Wbay
            rectangle(ax, 'Position',[x y DUCK_L DUCK_W], 'Curvature',0.1, ...
                      'FaceColor',[1.00 0.90 0.70], 'EdgeColor',[0.7 0.4 0.2]);
        end
    end

    xlim(ax, [-0.2, Lbay+0.2]); ylim(ax, [-0.2, Wbay+0.6]);
    axis(ax,'equal'); box(ax,'on');
    xlabel(ax,'length (in)'); ylabel(ax,'width (in)');
end
