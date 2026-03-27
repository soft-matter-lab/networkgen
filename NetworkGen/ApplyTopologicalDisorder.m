function [Bonds, Nvec] = ApplyTopologicalDisorder(obj, Atoms, Bonds, Nvec)
% -------------------------------------------------------------------------
% ApplyTopologicalDisorder
%   Probabilistically delete bonds to introduce topological defects into
%   a lattice network while respecting a minimum-degree constraint.
%
%   Only meaningful for hex_lattice geometry; the guard is enforced by
%   AddHeterogeneities before this function is called.
%
%   Reads from:
%     obj.architecture.lattice_disorder_level   [0, 1]
%         Controls expected deletions per node smoothly:
%             expected_del_per_node = disorder_level * lattice_max_del_per_node
%         At 0:  no bonds deleted.
%         At 1:  expected deletions per node = lattice_max_del_per_node.
%
%     obj.architecture.lattice_max_del_per_node
%         Upper bound on expected bond deletions per node.
%         Set to 0 to disable topological disorder entirely.
%         Default in architecture.m: 1
%
%     obj.architecture.lattice_min_degree_keep
%         Bonds are never deleted if either endpoint would fall below
%         this degree.  Prevents isolated or dangling nodes.
%         Default in architecture.m: 5
%
%   The per-bond deletion probability is derived from the expected
%   deletions per node and the current average degree:
%
%       lambda  = disorder_level * lattice_max_del_per_node
%       d_avg   = 2 * Nbonds / Natoms
%       p_del   = lambda / d_avg
%
%   Bonds are visited in random order.  A bond is deleted when:
%       rand < p_del  AND  deg(i) > min_degree_keep
%                     AND  deg(j) > min_degree_keep
%
%   Nvec is kept in sync with Bonds so downstream functions see a
%   consistent [Nbonds x ...] array.
%
% INPUT
%   obj   : network object
%   Atoms : [N x (5+MaxNbr)]  atom array  (read-only; not modified here)
%   Bonds : [M x 5]           bond array  [id | i | j | L0 | type]
%   Nvec  : [M x ...]         per-bond quantity array
%
% OUTPUT
%   Bonds : bond array with deleted bonds removed and IDs renumbered
%   Nvec  : filtered to match surviving bonds
% -------------------------------------------------------------------------

    disorder_level  = obj.architecture.lattice_disorder_level;
    max_del         = obj.architecture.lattice_max_del_per_node;
    min_keep        = obj.architecture.lattice_min_degree_keep;

    Natoms = size(Atoms, 1);
    Nbonds = size(Bonds, 1);

    if Nbonds == 0 || max_del <= 0 || disorder_level <= 0
        return;
    end

    % ------------------------------------------------------------------
    % Compute per-bond deletion probability from desired expected rate
    % ------------------------------------------------------------------
    lambda = disorder_level * max_del;   % expected deletions per node

    d_avg = 2.0 * Nbonds / max(Natoms, 1);
    if d_avg > 0
        p_del = lambda / d_avg;
    else
        p_del = 0;
    end

    % Hard clamp to avoid pathological behaviour at very high disorder
    p_del = max(0, min(0.9, p_del));

    if p_del <= 0
        return;
    end

    % ------------------------------------------------------------------
    % Build live degree vector from Bonds (independent of Atoms(:,5)
    % which may not be up to date after earlier heterogeneity steps)
    % ------------------------------------------------------------------
    deg = zeros(Natoms, 1);
    for k = 1:Nbonds
        ii = Bonds(k, 2);
        jj = Bonds(k, 3);
        if ii >= 1 && ii <= Natoms,  deg(ii) = deg(ii) + 1;  end
        if jj >= 1 && jj <= Natoms,  deg(jj) = deg(jj) + 1;  end
    end

    % ------------------------------------------------------------------
    % Visit bonds in random order and delete probabilistically
    % ------------------------------------------------------------------
    perm        = randperm(Nbonds);
    delete_flag = false(Nbonds, 1);
    n_deleted   = 0;

    for idx = 1:Nbonds
        k  = perm(idx);
        ii = Bonds(k, 2);
        jj = Bonds(k, 3);

        % Degree guard: never drop a node below min_keep
        if deg(ii) <= min_keep || deg(jj) <= min_keep
            continue;
        end

        if rand < p_del
            delete_flag(k) = true;
            deg(ii) = deg(ii) - 1;
            deg(jj) = deg(jj) - 1;
            n_deleted = n_deleted + 1;
        end
    end

    % ------------------------------------------------------------------
    % Apply deletion mask and renumber
    % ------------------------------------------------------------------
    keep_mask = ~delete_flag;
    Bonds     = Bonds(keep_mask, :);

    if ~isempty(Bonds)
        Bonds(:, 1) = (1:size(Bonds, 1)).';
    end

    % Keep Nvec in sync
    if ~isempty(Nvec)
        Nvec = Nvec(keep_mask, :);
    end

    obj.log.print(['   [ApplyTopologicalDisorder] Deleted %d / %d bonds ' ...
                   '(disorder_level=%.2f, p_del=%.4f, min_degree_keep=%d)\n'], ...
                  n_deleted, Nbonds, disorder_level, p_del, min_keep);

end
