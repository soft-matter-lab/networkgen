function [Atoms, Bonds, Nvec] = AddDefects(obj, Atoms, Bonds, Nvec)
% -------------------------------------------------------------------------
% AddDefects
%   Insert heterogeneous void defects into a polymer network.
%
%   Reads all settings from obj.defect, obj.domain, and obj.flags.
%   Void shape, placement, size distribution, and sparse/foam mode are
%   all controlled by obj.defect fields (see network.m for defaults).
%
%   NOTE: This function ONLY removes atoms and bonds that fall inside
%   voids.  It does NOT prune isolated atoms, remove disconnected
%   clusters, or renumber IDs.  All cleanup is deferred to CleanupNetwork.
%
% INPUT
%   obj   : network object (reads obj.defect, obj.domain, obj.flags)
%   Atoms : [N x (5+MaxNbr)]  atom array
%   Bonds : [M x 5]           bond array  [id | i | j | L0 | type]
%   Nvec  : [M x ...]         per-bond quantity array (e.g. Kuhn-N values)
%
% OUTPUT
%   Atoms : unchanged (atom positions / IDs not modified here)
%   Bonds : bonds surviving after void removal  (not yet renumbered)
%   Nvec  : Nvec rows filtered to match surviving bonds
% -------------------------------------------------------------------------

    %% ------------------------------------------------------------------
    %  0.  Early exit if defects are disabled
    %% ------------------------------------------------------------------
    if ~obj.flags.idefect
        fprintf('   [AddDefects] Skipped (flags.idefect = false)\n');
        return;
    end

    %% ------------------------------------------------------------------
    %  1.  Unpack settings from obj
    %% ------------------------------------------------------------------
    d = obj.defect;   % shorthand

    xlo = obj.domain.xlo;  xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;  yhi = obj.domain.yhi;
    Lx  = xhi - xlo;
    Ly  = yhi - ylo;
    domain_area = Lx * Ly;

    r_mean    = d.radius_mean;
    r_std     = d.radius_std;
    r_min     = d.radius_min;
    r_max     = d.radius_max;
    roughness = d.shape_roughness;
    n_modes   = max(1, round(d.shape_n_modes));
    wall_t    = d.wall_thickness;
    do_sparse = logical(d.sparse_network);

    %% ------------------------------------------------------------------
    %  2.  Determine number of voids
    %% ------------------------------------------------------------------
    if strcmpi(d.density_mode, 'area_frac')
        mean_void_area = pi * r_mean^2;
        n_voids = max(1, round(d.void_area_frac * domain_area / mean_void_area));
        fprintf('   [AddDefects] area_frac mode: targeting %.1f%% coverage -> %d voids\n', ...
            d.void_area_frac * 100, n_voids);
    else
        n_voids = max(1, round(d.n_voids));
        fprintf('   [AddDefects] count mode: %d voids requested\n', n_voids);
    end

    %% ------------------------------------------------------------------
    %  3.  Sample void radii
    %% ------------------------------------------------------------------
    radii = local_sample_radii(n_voids, d.size_dist, r_mean, r_std, r_min, r_max);

    %% ------------------------------------------------------------------
    %  4.  Compute placement zone (margin + clamp bands)
    %% ------------------------------------------------------------------
    margin  = d.margin_frac * r_mean;
    x_lo    = xlo + margin;   x_hi = xhi - margin;
    y_lo    = ylo + margin;   y_hi = yhi - margin;

    % Safety: if margin collapsed the range, revert to full domain
    if x_hi <= x_lo,  x_lo = xlo;  x_hi = xhi;  end
    if y_hi <= y_lo,  y_lo = ylo;  y_hi = yhi;  end

    % Clamp bands (defect-free strips at top and bottom)
    clamp_h   = d.clamp_thickness * Ly;
    clamp_ylo = ylo + clamp_h;
    clamp_yhi = yhi - clamp_h;

    if clamp_yhi <= clamp_ylo
        warning('AddDefects: clamp_thickness=%.2g leaves no active zone; clamps ignored.', ...
            d.clamp_thickness);
        clamp_h   = 0;
        clamp_ylo = ylo;
        clamp_yhi = yhi;
    end

    % Restrict void placement to the tighter of margin and clamp
    y_lo = max(y_lo, clamp_ylo);
    y_hi = min(y_hi, clamp_yhi);
    if y_hi <= y_lo
        y_lo = clamp_ylo;
        y_hi = clamp_yhi;
    end

    eff_Lx = x_hi - x_lo;
    eff_Ly = y_hi - y_lo;

    %% ------------------------------------------------------------------
    %  5.  Place void centers
    %% ------------------------------------------------------------------
    switch lower(d.center_distribution)

        case 'random'
            cand_x = x_lo + eff_Lx * rand(n_voids, 1);
            cand_y = y_lo + eff_Ly * rand(n_voids, 1);

        case 'uniform'
            % Jittered grid: one cell per void
            nx_g   = max(1, round(sqrt(n_voids * eff_Lx / eff_Ly)));
            ny_g   = max(1, ceil(n_voids / nx_g));
            cell_w = eff_Lx / nx_g;
            cell_h = eff_Ly / ny_g;

            [col_idx, row_idx] = meshgrid(0:nx_g-1, 0:ny_g-1);
            col_flat = col_idx(:);
            row_flat = row_idx(:);
            perm     = randperm(numel(col_flat));
            col_flat = col_flat(perm);
            row_flat = row_flat(perm);

            use    = min(n_voids, numel(col_flat));
            cand_x = zeros(n_voids, 1);
            cand_y = zeros(n_voids, 1);
            for v = 1:use
                cand_x(v) = x_lo + (col_flat(v) + rand) * cell_w;
                cand_y(v) = y_lo + (row_flat(v) + rand) * cell_h;
            end
            if use < n_voids
                cand_x(use+1:end) = x_lo + eff_Lx * rand(n_voids - use, 1);
                cand_y(use+1:end) = y_lo + eff_Ly * rand(n_voids - use, 1);
            end

        case 'clustered'
            % Thomas cluster process
            n_par  = d.n_cluster_parents;
            spread = d.cluster_spread;

            par_x  = x_lo + eff_Lx * rand(n_par, 1);
            par_y  = y_lo + eff_Ly * rand(n_par, 1);

            cand_x = zeros(n_voids, 1);
            cand_y = zeros(n_voids, 1);
            for v = 1:n_voids
                p = randi(n_par);
                for attempt = 1:50
                    cx_try = par_x(p) + spread * randn;
                    cy_try = par_y(p) + spread * randn;
                    if cx_try >= x_lo && cx_try <= x_hi && ...
                       cy_try >= y_lo && cy_try <= y_hi
                        break;
                    end
                end
                cand_x(v) = max(x_lo, min(x_hi, cx_try));
                cand_y(v) = max(y_lo, min(y_hi, cy_try));
            end

        otherwise
            warning('AddDefects: unknown center_distribution "%s"; using random.', ...
                d.center_distribution);
            cand_x = x_lo + eff_Lx * rand(n_voids, 1);
            cand_y = y_lo + eff_Ly * rand(n_voids, 1);
    end

    fprintf('   [AddDefects] center_distribution = %s\n', d.center_distribution);

    %% ------------------------------------------------------------------
    %  6.  Enforce no-overlap / bridge constraint (if requested)
    %% ------------------------------------------------------------------
    centers   = zeros(n_voids, 2);
    eff_radii = radii * (1 + roughness);   % worst-case boundary per void
    placed    = 0;

    % Accept logical, numeric, or string value for void_overlap
    vo = d.void_overlap;
    if ischar(vo) || isstring(vo)
        vo = ~(strcmpi(vo,'false') || strcmp(vo,'0'));
    end
    vo = logical(vo);

    if vo
        % Overlap allowed: accept all candidates
        centers = [cand_x, cand_y];
        placed  = n_voids;

    else
        % Overlap forbidden: enforce minimum edge-to-edge gap
        bridge           = max(0, d.bridge_width);
        max_tries_per_void = 200;

        for v = 1:n_voids
            accepted = false;
            for attempt = 1:max_tries_per_void
                if attempt == 1
                    cx_try = cand_x(v);
                    cy_try = cand_y(v);
                else
                    cx_try = x_lo + eff_Lx * rand;
                    cy_try = y_lo + eff_Ly * rand;
                end

                ok = true;
                for prev = 1:placed
                    d_edge = sqrt((cx_try - centers(prev,1))^2 + ...
                                  (cy_try - centers(prev,2))^2) ...
                             - eff_radii(v) - eff_radii(prev);
                    if d_edge < bridge
                        ok = false;
                        break;
                    end
                end

                if ok
                    placed = placed + 1;
                    centers(placed,1) = cx_try;
                    centers(placed,2) = cy_try;
                    accepted = true;
                    break;
                end
            end

            if ~accepted
                warning(['AddDefects: could not place void %d with ' ...
                         'bridge_width=%.2g after %d attempts; ' ...
                         'stopping at %d placed voids.'], ...
                         v, bridge, max_tries_per_void, placed);
                break;
            end
        end

        n_voids   = placed;
        radii     = radii(1:placed);
        eff_radii = eff_radii(1:placed);
        centers   = centers(1:placed, :);
    end

    fprintf('   [AddDefects] Placed %d voids  (void_overlap=%d)\n', n_voids, vo);

    %% ------------------------------------------------------------------
    %  7.  Generate per-void boundary roughness (polar Fourier perturbation)
    %
    %  Void boundary:  r(theta) = R * [1 + A * sum_k( w_k*cos(k*theta+phi_k) )]
    %  where A = shape_roughness, mode indices k = 2..n_modes+1 (k=1 skipped
    %  as it would merely translate the apparent center).
    %% ------------------------------------------------------------------
    mode_weights = zeros(n_voids, n_modes);
    mode_phases  = zeros(n_voids, n_modes);
    k_vals       = 2:(n_modes + 1);

    for v = 1:n_voids
        raw_amps           = (1 ./ k_vals) .* (0.5 + rand(1, n_modes));
        mode_weights(v,:)  = raw_amps / sum(raw_amps);
        mode_phases(v,:)   = 2 * pi * rand(1, n_modes);
    end

    %% ------------------------------------------------------------------
    %  8.  Classify atoms as inside void / near boundary
    %% ------------------------------------------------------------------
    natom        = size(Atoms, 1);
    in_void       = false(natom, 1);
    near_boundary = false(natom, 1);   % only meaningful in sparse mode

    atom_x = Atoms(:, 2);
    atom_y = Atoms(:, 3);

    for v = 1:n_voids
        cx_v = centers(v, 1);
        cy_v = centers(v, 2);
        R_v  = radii(v);
        wts  = mode_weights(v, :);
        phis = mode_phases(v, :);

        dx       = atom_x - cx_v;
        dy       = atom_y - cy_v;
        dist2    = dx .* dx + dy .* dy;

        % Broad-phase: only examine atoms within worst-case outer radius
        outer_R   = (1.0 + roughness) * R_v + wall_t;
        candidate = find(dist2 < outer_R * outer_R);
        if isempty(candidate),  continue;  end

        % Compute local (rough) void radius at each candidate's angle
        theta_c  = atan2(dy(candidate), dx(candidate));

        k_mat    = repmat(k_vals,     numel(candidate), 1);
        phi_mat  = repmat(phis,       numel(candidate), 1);
        w_mat    = repmat(wts,        numel(candidate), 1);
        th_mat   = repmat(theta_c,    1, n_modes);

        perturb  = sum(w_mat .* cos(k_mat .* th_mat + phi_mat), 2);
        r_local  = R_v * max(0.1, 1 + roughness * perturb);

        dist_c   = sqrt(dist2(candidate));
        inside   = dist_c < r_local;

        in_void(candidate(inside)) = true;

        if do_sparse
            near = ~inside & (dist_c < r_local + wall_t);
            near_boundary(candidate(near)) = true;
        end
    end

    % Protect clamp-band atoms — they must never be removed
    if clamp_h > 0
        in_clamp               = (atom_y <= clamp_ylo) | (atom_y >= clamp_yhi);
        in_void(in_clamp)       = false;
        near_boundary(in_clamp) = true;
    end

    % Final removal mask
    if do_sparse
        remove_mask = in_void | ~near_boundary;
        fprintf('   [AddDefects] sparse_network mode  (wall_thickness=%.2g)\n', wall_t);
    else
        remove_mask = in_void;
    end

    n_removed = sum(remove_mask);
    fprintf('   [AddDefects] Marking %d / %d atoms for removal (%.1f%%)\n', ...
        n_removed, natom, 100.0 * n_removed / max(natom, 1));

    %% ------------------------------------------------------------------
    %  9.  Remove atoms inside voids
    %% ------------------------------------------------------------------
    Atoms = Atoms(~remove_mask, :);

    if isempty(Atoms)
        warning(['AddDefects: ALL atoms were removed. ' ...
                 'Void radii or count may be too large for this domain.']);
        Bonds = zeros(0, size(Bonds, 2));
        Nvec  = [];
        return;
    end

    surviving_ids = Atoms(:, 1);

    %% ------------------------------------------------------------------
    %  10.  Remove bonds whose endpoints were removed
    %% ------------------------------------------------------------------
    if isempty(Bonds)
        Bonds = zeros(0, size(Bonds, 2));
        Nvec  = [];
        fprintf('   [AddDefects] No bonds to filter\n');
        return;
    end

    keep1    = ismember(Bonds(:, 2), surviving_ids);
    keep2    = ismember(Bonds(:, 3), surviving_ids);
    keep     = keep1 & keep2;

    n_bonds_rm = sum(~keep);
    Bonds = Bonds(keep, :);
    fprintf('   [AddDefects] Removed %d / %d bonds incident to void atoms\n', ...
        n_bonds_rm, n_bonds_rm + size(Bonds, 1));

    %% ------------------------------------------------------------------
    %  11.  Filter Nvec to match surviving bonds
    %% ------------------------------------------------------------------
    if ~isempty(Nvec)
        Nvec = Nvec(keep, :);
    end

    fprintf(['   [AddDefects] Done. ' ...
             '%d atoms remain (-%d), %d bonds remain (-%d).\n'], ...
        size(Atoms, 1),  n_removed, ...
        size(Bonds, 1),  n_bonds_rm);

end   % end main function


%% =======================================================================
%  LOCAL HELPER:  sample_radii
%% =======================================================================
function radii = local_sample_radii(n, size_dist, r_mean, r_std, r_min, r_max)
% Draw n void radii from the requested distribution.
%   'fixed'       All radii equal r_mean
%   'gaussian'    Normal(r_mean, r_std), clipped to [r_min, r_max]
%   'exponential' Exponential with mean=r_mean, shifted to start at r_min,
%                 clipped at r_max.  No Statistics Toolbox required.

    switch lower(size_dist)

        case 'fixed'
            radii = r_mean * ones(n, 1);

        case 'gaussian'
            radii = r_mean + r_std * randn(n, 1);
            radii = max(r_min, min(r_max, radii));

        case 'exponential'
            % Inverse-transform sampling: X = -mu*log(U) ~ Exp(mu)
            adjusted_mean = max(r_mean - r_min, 1e-10);
            U     = rand(n, 1);
            U     = max(U, 1e-15);
            radii = r_min + (-adjusted_mean * log(U));
            radii = min(radii, r_max);

        otherwise
            warning('AddDefects:sample_radii:unknownDist', ...
                'Unknown size_dist "%s"; defaulting to fixed.', size_dist);
            radii = r_mean * ones(n, 1);
    end
end
