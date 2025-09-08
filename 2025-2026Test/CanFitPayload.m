function [fits, info] = CanFitPayload(plane, p, c, payloadBayDims)
% Uses a simple length-partition packing:
%   [puck zone along length] + [duck zone along length] <= bay length
% INPUTS
%   plane            : AirplaneClass/struct (uses fuselage dims if dims not given)
%   p                : # ducks (passengers)
%   c                : # pucks (cargo)
%   payloadBayDims   : (optional) struct with fields (inches):
%                     .length_in, .width_in, .height_
% OUTPUTS
%   fits   : true/false
%   info   : struct with fields describing capacities and required lengths

    % ---- Duck & puck "bounding box" dimensions (inches) -------------------
    DUCK_W = 2.3; DUCK_L = 2.5; DUCK_H = 2.5;   % upright, single compartment
    PUCK_D = 3.0; PUCK_T = 1.0;                 % flat puck (thickness vertical)

    % ---- Bay dimensions (inches) ------------------------------------------
    if nargin >= 4 && ~isempty(payloadBayDims) ...
            && all(isfield(payloadBayDims, {'length_in','width_in','height_in'}))
        L_in = payloadBayDims.length_in;
        W_in = payloadBayDims.width_in;
        H_in = payloadBayDims.height_in;
    else
        % Derive from fuselage, with conservative interior margins
        wallClr_in = 0.25;  % each side
        L_in = max(10, 0.40 * (plane.fuselage.length / 0.0254));   % use ~40% fuselage as payload bay
        W_in = max(3.0, (plane.fuselage.width  / 0.0254) - 2*wallClr_in);
        H_in = max(3.0, (plane.fuselage.height / 0.0254) - 2*wallClr_in);
    end

    % sanity clamps
    L_in = max(0.1, L_in); W_in = max(0.1, W_in); H_in = max(0.1, H_in);

    % ---- Discrete capacities per "slice" along length ---------------------
    % Ducks (upright): footprint (DUCK_W x DUCK_L), height DUCK_H
    cap_duck_per_slice = floor(W_in / DUCK_W) * floor(H_in / DUCK_H);  % # ducks per 2.5" length slice
    if cap_duck_per_slice < 1, cap_duck_per_slice = 0; end

    % Pucks (lying flat): footprint (PUCK_D x PUCK_D), height PUCK_T
    cap_puck_per_slice = floor(W_in / PUCK_D) * floor(H_in / PUCK_T);  % # pucks per 3.0" length slice
    if cap_puck_per_slice < 1, cap_puck_per_slice = 0; end

    % ---- Required lengths for requested quantities ------------------------
    if cap_duck_per_slice == 0 && p > 0
        L_duck_req = inf;
    else
        L_duck_req = ceil(p / max(cap_duck_per_slice,1)) * DUCK_L;
    end

    if cap_puck_per_slice == 0 && c > 0
        L_puck_req = inf;
    else
        L_puck_req = ceil(c / max(cap_puck_per_slice,1)) * PUCK_D;
    end

    L_total_req = L_duck_req + L_puck_req;
    fits = (L_total_req <= L_in);

    % ---- Pack results ------------------------------------------------------
    info = struct();
    info.bay.length_in = L_in; info.bay.width_in = W_in; info.bay.height_in = H_in;
    info.duck.cap_per_slice = cap_duck_per_slice; info.duck.slice_len_in = DUCK_L; info.duck.len_req_in = L_duck_req;
    info.puck.cap_per_slice = cap_puck_per_slice; info.puck.slice_len_in = PUCK_D; info.puck.len_req_in = L_puck_req;
    info.total_required_len_in = L_total_req;
end
