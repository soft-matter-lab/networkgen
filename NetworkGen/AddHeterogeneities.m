function [Atoms, Bonds, Nvec] = AddHeterogeneities(obj, Atoms, Bonds, Nvec)
% -------------------------------------------------------------------------
% AddHeterogeneities
%   Top-level dispatcher for all network heterogeneity operations.
%   Drop-in replacement for the bare AddDefects call in generateNetwork.
%
%   Runs up to three sub-operations in order:
%
%     1. Void defects  (any geometry)
%            Controlled by obj.flags.idefect and obj.defect.*
%            Removes atoms/bonds that fall inside void regions.
%            -> calls AddDefects
%
%     2. Geometric lattice disorder  (hex_lattice only)
%            Controlled by obj.architecture.lattice_disorder_level > 0
%            Randomly perturbs atom positions and recomputes L0.
%            -> calls ApplyGeometricDisorder
%
%     3. Topological lattice disorder  (hex_lattice only)
%            Controlled by obj.architecture.lattice_disorder_level > 0
%                      AND obj.architecture.lattice_max_del_per_node > 0
%            Probabilistically deletes bonds respecting min_degree_keep.
%            -> calls ApplyTopologicalDisorder
%
%   Geometric disorder is applied before topological disorder so that
%   realistic L0 lengths are in place before any bonds are deleted.
%
% INPUT
%   obj   : network object
%   Atoms : [N x (5+MaxNbr)]  atom array
%   Bonds : [M x 5]           bond array  [id | i | j | L0 | type]
%   Nvec  : [M x ...]         per-bond quantity array
%
% OUTPUT
%   Atoms : modified atom array  (positions and/or count changed)
%   Bonds : modified bond array  (count and/or L0 changed)
%   Nvec  : filtered to match surviving bonds
% -------------------------------------------------------------------------

    isHex = strcmpi(obj.architecture.geometry, 'hex_lattice');

    % ------------------------------------------------------------------
    % 1.  Void defects  (geometry-agnostic)
    % ------------------------------------------------------------------
    [Atoms, Bonds, Nvec] = AddDefects(obj, Atoms, Bonds, Nvec);

    % ------------------------------------------------------------------
    % 2.  Geometric disorder  (hex_lattice only)
    % ------------------------------------------------------------------
    if isHex && obj.architecture.lattice_disorder_level > 0
        [Atoms, Bonds] = ApplyGeometricDisorder(obj, Atoms, Bonds);
    else
        if ~isHex && obj.architecture.lattice_disorder_level > 0
            obj.log.print(['   [AddHeterogeneities] Geometric disorder skipped ' ...
                           '(only supported for hex_lattice geometry)\n']);
        end
    end

    % ------------------------------------------------------------------
    % 3.  Topological disorder  (hex_lattice only)
    % ------------------------------------------------------------------
    if isHex && obj.architecture.lattice_disorder_level > 0 && ...
                obj.architecture.lattice_max_del_per_node > 0
        [Bonds, Nvec] = ApplyTopologicalDisorder(obj, Atoms, Bonds, Nvec);
    end

end
