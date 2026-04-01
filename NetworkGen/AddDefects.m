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
        obj.log.print('   [AddDefects] Skipped (flags.idefect = false)\n');
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
        obj.log.print('   [AddDefects] area_frac mode: targeting %.1f%% coverage -> %d voids\n', ...
            d.void_area_frac * 100, n_voids);
    else
        n_voids = max(1, round(d.n_voids));
        obj.log.print('   [AddDefects] count mode: %d voids requested\n', n_voids);
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

    obj.log.print('   [AddDefects] center_distribution = %s\n', d.center_distribution);

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

    obj.log.print('   [AddDefects] Placed %d voids  (void_overlap=%d)\n', n_voids, vo);

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

    % Protect clamp-band atoms (y-boundaries) — they must never be removed
    if clamp_h > 0
        in_clamp               = (atom_y <= clamp_ylo) | (atom_y >= clamp_yhi);
        in_void(in_clamp)       = false;
        near_boundary(in_clamp) = true;
    end

    % Protect x-boundary atoms in sparse mode.
    % When margin_frac is large, no voids are placed near the x-edges so
    % those atoms would never be marked near_boundary and get incorrectly removed.
    if do_sparse
        in_x_wall                = (atom_x <= xlo + wall_t) | (atom_x >= xhi - wall_t);
        in_void(in_x_wall)       = false;
        near_boundary(in_x_wall) = true;
    end

    % Final removal mask
    if do_sparse
        remove_mask = in_void | ~near_boundary;
        obj.log.print('   [AddDefects] sparse_network mode  (wall_thickness=%.2g)\n', wall_t);
    else
        remove_mask = in_void;
    end

    n_removed = sum(remove_mask);
    obj.log.print('   [AddDefects] Marking %d / %d atoms for removal (%.1f%%)\n', ...
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
        obj.log.print('   [AddDefects] No bonds to filter\n');
        return;
    end

    keep1    = ismember(Bonds(:, 2), surviving_ids);
    keep2    = ismember(Bonds(:, 3), surviving_ids);
    keep     = keep1 & keep2;

    n_bonds_rm = sum(~keep);
    Bonds = Bonds(keep, :);
    obj.log.print('   [AddDefects] Removed %d / %d bonds incident to void atoms\n', ...
        n_bonds_rm, n_bonds_rm + size(Bonds, 1));

    %% ------------------------------------------------------------------
    %  11.  Filter Nvec to match surviving bonds
    %% ------------------------------------------------------------------
    if ~isempty(Nvec)
        Nvec = Nvec(keep, :);
    end

    obj.log.print(['   [AddDefects] Done. ' ...
             '%d atoms remain (-%d), %d bonds remain (-%d).\n'], ...
        size(Atoms, 1),  n_removed, ...
        size(Bonds, 1),  n_bonds_rm);

    %% ------------------------------------------------------------------
    %  12.  Secondary thinning pass
    %       Reduces density in over-dense regions (bridges, untouched patches)
    %       by probabilistic atom removal weighted by local density.
    %       Articulation point check (optional) prevents disconnection.
    %% ------------------------------------------------------------------
    if isfield(d, 'thinning') && logical(d.thinning) && ~isempty(Atoms) && ~isempty(Bonds)
        [Atoms, Bonds, Nvec] = local_thin_network(obj, Atoms, Bonds, Nvec);
    end

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


%% =======================================================================
%  LOCAL HELPER:  local_thin_network
%
%  Secondary thinning pass after primary void placement.
%  Probabilistically removes atoms in over-dense regions — bridges between
%  voids and untouched dense patches — while optionally protecting
%  articulation points (cut vertices) to avoid disconnecting the network.
%
%  Reads from obj.defect:
%    thinning_radius      - neighbourhood radius for local density estimate
%                           (default: 2.5 × expected node spacing)
%    thinning_target_frac - target fraction of atoms to keep in dense regions
%                           (default: 0.4)
%    thinning_min_keep    - floor: minimum keep fraction anywhere
%                           (default: 0.1)
%    thinning_protect_art - if true, never remove articulation points
%                           (default: true)
%    thinning_recompute_k - recompute articulation points every k removals
%                           (default: 20)
%% =======================================================================
function [Atoms, Bonds, Nvec] = local_thin_network(obj, Atoms, Bonds, Nvec)

    d = obj.defect;

    % ── Parameters with defaults ────────────────────────────────────────
    rho_atom   = obj.architecture.rho_atom;
    avg_spacing = 1.0 / sqrt(max(rho_atom, 1e-10));

    if isfield(d, 'thinning_radius')
        r_thin = d.thinning_radius;
    else
        r_thin = 2.5 * avg_spacing;
    end

    if isfield(d, 'thinning_target_frac')
        target_frac = d.thinning_target_frac;
    else
        target_frac = 0.4;
    end

    if isfield(d, 'thinning_min_keep')
        min_keep = d.thinning_min_keep;
    else
        min_keep = 0.1;
    end

    if isfield(d, 'thinning_protect_art')
        protect_art = logical(d.thinning_protect_art);
    else
        protect_art = true;
    end

    if isfield(d, 'thinning_recompute_k')
        recompute_k = max(1, round(d.thinning_recompute_k));
    else
        recompute_k = 20;
    end

    obj.log.print('   [AddDefects:Thin] Starting secondary thinning pass\n');
    obj.log.print('   [AddDefects:Thin] r_thin=%.1f  target_frac=%.2f  protect_art=%d\n', ...
        r_thin, target_frac, protect_art);

    natom_before = size(Atoms, 1);

    % ── Atom positions ───────────────────────────────────────────────────
    atom_ids = Atoms(:, 1);
    ax       = Atoms(:, 2);
    ay       = Atoms(:, 3);
    natom    = numel(atom_ids);

    % Build a fast index: map atom ID -> row index
    max_id      = max(atom_ids);
    id_to_row   = zeros(max_id, 1, 'int32');
    for i = 1:natom
        id_to_row(atom_ids(i)) = i;
    end

    % ── Expected atoms in a circle of radius r_thin ──────────────────────
    expected_in_r = rho_atom * pi * r_thin^2;
    expected_in_r = max(expected_in_r, 1.0);   % avoid divide-by-zero

    % ── Compute local density for each atom ─────────────────────────────
    local_density = zeros(natom, 1);
    for i = 1:natom
        dx2      = (ax - ax(i)).^2 + (ay - ay(i)).^2;
        n_in_r   = sum(dx2 <= r_thin^2) - 1;   % exclude self
        local_density(i) = n_in_r / expected_in_r;
    end

    % ── Compute removal probability per atom ────────────────────────────
    % p_remove = 0                     when density <= 1 (at or below average)
    % p_remove = 1 - target_frac/rho   when density > 1  (over-dense)
    % clamped so keep fraction never drops below min_keep
    p_remove = zeros(natom, 1);
    over_dense = local_density > 1.0;
    p_remove(over_dense) = 1.0 - target_frac ./ local_density(over_dense);
    p_remove = min(p_remove, 1.0 - min_keep);
    p_remove = max(p_remove, 0.0);

    % Sort atoms: remove highest-density first so articulation point
    % checks reflect the most-congested regions first
    [~, sort_order] = sort(local_density, 'descend');

    % ── Initial articulation point computation ───────────────────────────
    alive = true(natom, 1);   % tracks which atoms are still present
    bond_i = Bonds(:, 2);
    bond_j = Bonds(:, 3);

    if protect_art
        is_art = local_find_articulation_points(atom_ids, bond_i, bond_j, alive, id_to_row);
    else
        is_art = false(natom, 1);
    end

    % ── Main removal loop ────────────────────────────────────────────────
    n_removed   = 0;
    since_last  = 0;

    for k = 1:natom
        i = sort_order(k);

        if ~alive(i),        continue; end
        if p_remove(i) <= 0, continue; end
        if protect_art && is_art(i), continue; end

        % Probabilistic removal
        if rand() < p_remove(i)
            alive(i)   = true;   % mark for removal
            alive(i)   = false;
            n_removed  = n_removed + 1;
            since_last = since_last + 1;

            % Recompute articulation points every recompute_k removals
            if protect_art && (since_last >= recompute_k)
                is_art     = local_find_articulation_points(atom_ids, bond_i, bond_j, alive, id_to_row);
                since_last = 0;
            end
        end
    end

    % Final articulation recompute to be safe before applying
    if protect_art && since_last > 0
        is_art = local_find_articulation_points(atom_ids, bond_i, bond_j, alive, id_to_row);
        % Un-remove any articulation points that slipped through
        for k = 1:natom
            if ~alive(k) && is_art(k)
                alive(k)  = true;
                n_removed = n_removed - 1;
            end
        end
    end

    % ── Apply removal ────────────────────────────────────────────────────
    surviving_ids = atom_ids(alive);
    Atoms = Atoms(alive, :);

    keep_bond = ismember(bond_i, surviving_ids) & ismember(bond_j, surviving_ids);
    n_bonds_rm = sum(~keep_bond);
    Bonds = Bonds(keep_bond, :);

    if ~isempty(Nvec)
        Nvec = Nvec(keep_bond, :);
    end

    obj.log.print(['   [AddDefects:Thin] Removed %d / %d atoms  ' ...
                   '(%d bonds removed).  %d atoms remain.\n'], ...
        n_removed, natom_before, n_bonds_rm, size(Atoms, 1));

end


%% =======================================================================
%  LOCAL HELPER:  local_find_articulation_points
%
%  Finds all articulation points (cut vertices) in the graph defined by
%  atom_ids / bond endpoints using Tarjan's iterative DFS algorithm.
%  Runs in O(N + M).
%
%  Inputs
%    atom_ids : [N x 1]  atom ID values
%    bond_i   : [M x 1]  bond endpoint atom IDs (column 2 of Bonds)
%    bond_j   : [M x 1]  bond endpoint atom IDs (column 3 of Bonds)
%    alive    : [N x 1]  logical — which atoms are currently present
%    id_to_row: lookup array, id_to_row(atom_id) = row index in atom_ids
%
%  Output
%    is_art   : [N x 1] logical — true if atom is an articulation point
%% =======================================================================
function is_art = local_find_articulation_points(atom_ids, bond_i, bond_j, alive, id_to_row)

    natom  = numel(atom_ids);
    is_art = false(natom, 1);

    % Build adjacency list for alive atoms only
    % Only include bonds where both endpoints are alive
    alive_ids = atom_ids(alive);
    alive_set = false(numel(id_to_row), 1);
    alive_set(alive_ids) = true;

    valid_bond = alive_set(bond_i) & alive_set(bond_j);
    vi = bond_i(valid_bond);
    vj = bond_j(valid_bond);

    % Map to row indices
    ri = id_to_row(vi);
    rj = id_to_row(vj);

    % Build adjacency list
    adj = cell(natom, 1);
    for k = 1:numel(ri)
        a = ri(k);  b = rj(k);
        if alive(a) && alive(b)
            adj{a}(end+1) = b;
            adj{b}(end+1) = a;
        end
    end

    % Tarjan's articulation point algorithm — iterative DFS
    disc    = zeros(natom, 1);   % discovery time
    low     = zeros(natom, 1);   % low value
    parent  = zeros(natom, 1);   % parent in DFS tree
    visited = false(natom, 1);
    timer   = 0;

    alive_rows = find(alive);

    for start_idx = 1:numel(alive_rows)
        s = alive_rows(start_idx);
        if visited(s), continue; end

        % Iterative DFS using explicit stack
        % Stack entries: [node, parent, child_iterator_index]
        stack      = zeros(natom, 3, 'int32');
        stack_top  = 1;
        stack(1,:) = [s, 0, 1];

        timer      = timer + 1;
        disc(s)    = timer;
        low(s)     = timer;
        visited(s) = true;
        parent(s)  = -1;   % root has no parent

        while stack_top > 0
            node      = stack(stack_top, 1);
            par       = stack(stack_top, 2);
            child_idx = stack(stack_top, 3);

            nbrs = adj{node};

            if child_idx > numel(nbrs)
                % Done with this node — pop and update parent's low
                stack_top = stack_top - 1;
                if par > 0
                    low(par) = min(low(par), low(node));

                    % Articulation point check (non-root)
                    if parent(par) > 0 && low(node) >= disc(par)
                        is_art(par) = true;
                    end
                end
            else
                % Advance child iterator
                stack(stack_top, 3) = child_idx + 1;
                nb = nbrs(child_idx);

                if ~visited(nb)
                    timer       = timer + 1;
                    disc(nb)    = timer;
                    low(nb)     = timer;
                    visited(nb) = true;
                    parent(nb)  = node;

                    stack_top = stack_top + 1;
                    stack(stack_top, :) = [nb, node, 1];

                elseif nb ~= par
                    low(node) = min(low(node), disc(nb));
                    % Update stack entry for back-edge
                    stack(stack_top, 2) = par;  % keep parent unchanged
                    % low update already done — re-store node entry
                    stack(stack_top, 1) = node;
                end
            end
        end

        % Root articulation point check: root is an art point if it has
        % more than one child in the DFS tree
        root_children = sum(parent(alive_rows) == s);
        if root_children > 1
            is_art(s) = true;
        end
    end

end