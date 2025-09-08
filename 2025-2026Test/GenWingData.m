function [wings, airfoilTable] = GenWingData(Aspect_Ratios, span_m, varargin)
% GenWingData
% Builds wing lookup tensor:
%   wings(airfoilIdx, param, ARidx, spanIdx)
% with columns:
%   1 C_l(3°), 2 C_lmax(15°), 3 C_di, 4 C_l_flaps, 5 weight(kg),
%   6 chord(m), 7 planformArea(m^2), 8 surfaceArea(m^2),
%   9 airfoil ID/1000 (legacy), 10 thickness(m)
% Outputs:
%   wings         — 4-D lookup as above
%   airfoilTable  — struct array with fields: name,id,alpha0_deg,t_over_c,xt_over_c
% Notes:
% - Lifting-line uses alpha0_deg, alphac=3°, alphat=15°.
% - cl3/cd3 kept only for reference.
% - HQ 3.0/12 cd reference corrected (~0.012).

AR  = Aspect_Ratios(:).';
span= span_m(:).';

% ----- Airfoil DB (extend/edit here) --------------------------------------
db = [ ...
%  name        id     alpha0  t/c    x_t/c  cl@3°   cd@3°
struct('name','SD7062',   'id',7062, 'alpha0_deg',-4.235,'t_over_c',0.1398,'xt_over_c',0.2715,'cl3',0.77485,'cd3',0.01249)
struct('name','SD8040',   'id',8040, 'alpha0_deg',-2.714,'t_over_c',0.1000,'xt_over_c',0.2933,'cl3',0.64702,'cd3',0.00997)
struct('name','Eppler66', 'id',66,   'alpha0_deg',-3.589,'t_over_c',0.1013,'xt_over_c',0.3145,'cl3',0.86236,'cd3',0.01106)
struct('name','Eppler374','id',374,  'alpha0_deg',-2.012,'t_over_c',0.1091,'xt_over_c',0.3434,'cl3',0.63729,'cd3',0.01185)
struct('name','Goettingen412','id',412,'alpha0_deg',-5.805,'t_over_c',0.1313,'xt_over_c',0.2973,'cl3',0.86640,'cd3',0.01340)
struct('name','Goettingen433','id',433,'alpha0_deg',-5.184,'t_over_c',0.1741,'xt_over_c',0.2913,'cl3',0.91441,'cd3',0.01548)
struct('name','NACA2414','id',2414, 'alpha0_deg',-2.123,'t_over_c',0.1400,'xt_over_c',0.2953,'cl3',0.64901,'cd3',0.01163)
struct('name','HQ3.0/12','id',3003, 'alpha0_deg',-3.591,'t_over_c',0.1198,'xt_over_c',0.3504,'cl3',0.71923,'cd3',0.01200)  % fixed cd
% ---- Added RC favorites ---------------------------------------------------
struct('name','S1223',    'id',1223, 'alpha0_deg',-2.8,  't_over_c',0.1220,'xt_over_c',0.33,  'cl3',0.95,  'cd3',0.0120)
struct('name','MH32',     'id',3200, 'alpha0_deg',-2.0,  't_over_c',0.1200,'xt_over_c',0.33,  'cl3',0.80,  'cd3',0.0105)
struct('name','Eppler205','id',205,  'alpha0_deg',-1.5,  't_over_c',0.1000,'xt_over_c',0.30,  'cl3',0.70,  'cd3',0.0110)
struct('name','S3021',    'id',3021, 'alpha0_deg',-2.5,  't_over_c',0.1250,'xt_over_c',0.32,  'cl3',0.80,  'cd3',0.0115)
struct('name','NACA23012','id',23012,'alpha0_deg',-2.0,  't_over_c',0.1200,'xt_over_c',0.30,  'cl3',0.70,  'cd3',0.0120)
];

% (The cl3/cd3 above are reasonable RC-scale references; replace with your
% XFOIL or wind-tunnel polars at your Reynolds for highest fidelity.)

% ----- Optional selection by names ----------------------------------------
if ~isempty(varargin)
    wantNames = string(varargin{1});
    keep = false(size(db));
    for k = 1:numel(db)
        keep(k) = any(strcmpi(db(k).name, wantNames));
    end
    if ~any(keep)
        error('GenWingData:NoMatch','Requested airfoil names not found.');
    end
    db = db(keep);
end

% ----- Physical/material constants (balsa + CF spars like your original) ---
rho_balsa = 160.3;         % kg/m^3 (aircraft balsa)
w_spar    = 0.1711;        % kg/m (CF spars from McMaster ref)
in2m      = 0.0254;

% ----- Lifting-line setup --------------------------------------------------
alphac = 3*pi/180;            % cruise
alphat = 15*pi/180;           % takeoff
theta  = [pi/4 pi/2 3*pi/4 pi];
n      = [1 3 5 7];

NAF = numel(db);
wings = zeros(NAF,10,numel(AR),numel(span));

for iAR = 1:numel(AR)
    for iaf = 1:NAF
        mu = pi/(2*AR(iAR));
        alpha0 = db(iaf).alpha0_deg*pi/180;
        alphae   = alphac - alpha0;
        alphamax = alphat - alpha0;

        % Precompute psi once per AR
        psi = zeros(4,4);
        for z=1:4
            s = sin(theta(z));
            psi(z,:) = [s*(mu+s), sin(3*theta(z))*(3*mu+s), sin(5*theta(z))*(5*mu+s), sin(7*theta(z))*(7*mu+s)];
        end

        % zeta vectors for cruise/takeoff
        zeta    = mu*alphae   * sin(theta(:));
        zetamax = mu*alphamax * sin(theta(:));

        A    = psi\zeta;
        Amax = psi\zetamax;

        C_l    = A(1)*pi*AR(iAR);
        C_lmax = Amax(1)*pi*AR(iAR);
        C_di   = pi*AR(iAR)*dot(n, A.^2);

        for ispan = 1:numel(span)
            chord = span(ispan)/AR(iAR);
            t     = chord*db(iaf).t_over_c;
            planA = chord*span(ispan);
            surfA = 2*planA;

            % Balsa weight estimate (your legacy formula, incl. ribs/stringers)
            ax      = (chord*db(iaf).xt_over_c*t/2) + (chord*(1-db(iaf).xt_over_c)*t/2);
            v_balsa = (span(ispan)*ax*1/32) + ((1/8*in2m)^2 * span(ispan)*8);
            Wkg     = (rho_balsa*v_balsa) + 2*w_spar*span(ispan);
            Wkg     = 1.30 * Wkg;   % epoxy/film fudge factor

            % Simple flap increment (same logic as your original)
            deltaClmax = (0.95)*(0.58)*(0.28)*1.15;
            C_lflaps   = C_lmax + deltaClmax;

            % Populate table
            wings(iaf,1,iAR,ispan)  = C_l;
            wings(iaf,2,iAR,ispan)  = C_lmax;
            wings(iaf,3,iAR,ispan)  = C_di;          % NOTE: stored as "cd" later
            wings(iaf,4,iAR,ispan)  = C_lflaps;
            wings(iaf,5,iAR,ispan)  = Wkg;
            wings(iaf,6,iAR,ispan)  = chord;
            wings(iaf,7,iAR,ispan)  = planA;
            wings(iaf,8,iAR,ispan)  = surfA;
            wings(iaf,9,iAR,ispan)  = db(iaf).id/1000;   % legacy name code
            wings(iaf,10,iAR,ispan) = t;
        end
    end
end

% Optional info table back to caller
if nargout > 1
    airfoilTable = rmfield(db, {'cl3','cd3'}); % keep essentials compact
end
end
