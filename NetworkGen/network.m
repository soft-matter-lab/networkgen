classdef network < handle

properties

    % Define properites and set default parameters

    %%% Flags
    flags = struct(...
        'isave', true, ...
        'iplot', true, ...
        'savemode', true, ...
        'imanualseed', false ...
    );

    %%% Domain
    domain = struct(...
        'b',  1.6, ...
        'Lx', 10, ...
        'Ly', 10, ...
        'Lz', 10, ...
        'scale', 1, ...
        'boundary', 'fixed', ...
        'xlo', [], ...
        'xhi', [], ...
        'ylo', [], ...
        'yhi', [], ...
        'zlo', [], ...
        'zhi', [], ...
        'Max_atom', [], ...
        'Max_bond', [], ...
        'node_scatter_max_tries', [], ...
        'bond_global_try_limit', [], ...
        'max_attempts_without_progress', [], ...
        'max_tries_per_node_sample', [] ...
    );

    %%% Architecture subclass
    arch architecture 

    %%% Per/atom
    peratom = struct(...
        'Max_peratom_bond',[] ...
    );

    %%% Per/bond subclass
    perbond bondstyle

    %%% Defect
    defect = struct(...
    );

    %%% Potential
    pot = struct(...
    );

    %%% log
    log = struct(...
    );

end

methods

    % Explicitly initialize subclass objects in the constructor
    function obj = network()
        obj.arch = architecture();
        obj.perbond = bondstyle();
    end

    function [] = generateNetwork(obj)

        %%% Loop over replicates
        for ii=1 %Nreps

            % ---------------------------------------------------------
            % 1. Prepare replicate-specific information
            % ---------------------------------------------------------
            obj = localPrepareReplicate(obj, ii, Nreplicates);

            % ---------------------------------------------------------
            % 2. Construct domain
            % ---------------------------------------------------------
            SetupDomain(obj);
            % New version of SetupDomain should read from obj.domain, % obj.arch, ...
            
            % ---------------------------------------------------------
            % 3. Add atoms
            % ---------------------------------------------------------
            Atoms = AddAtoms(obj);

            % ---------------------------------------------------------
            % 4. Assign per/atom
            % ---------------------------------------------------------
            Atoms = AssignPerAtom(obj, Atoms);
            % AssignPerAtom should read per-atom settings (if any) from obj
            % Decides internally to do random or hex based on geometry flag

            % ---------------------------------------------------------
            % 5. Add bonds
            % ---------------------------------------------------------
            [Atoms, Bonds] = AddBonds(obj, Atoms);

            % ---------------------------------------------------------
            % 6. Assign per/bond
            % ---------------------------------------------------------
            Nvec = AssignPerBond(obj, Bonds, Atoms);
            % AssignPerBond should read obj.perbond.* settings and assign

            % ---------------------------------------------------------
            % 7. Add defects
            % ---------------------------------------------------------
            [Atoms, Bonds] = AddDefects(obj, Atoms, Bonds, Nvec);
            % AddDefects should read obj.defect

            % ---------------------------------------------------------
            % 8. Clean-up network
            % ---------------------------------------------------------
            [Atoms, Bonds, Nvec] = CleanupNetwork(obj, Atoms, Bonds, Nvec);
            % CleanupNetwork can prune isolated nodes, rebuild connectivity,
            % update Nvec, etc.


            % ---------------------------------------------------------
            % 9. Construct local density potential
            % ---------------------------------------------------------
            LDpot = ConstructLDPotential(obj, Atoms, Bonds, Nvec);

            % ---------------------------------------------------------
            % 10. Scale domain if needed
            % ---------------------------------------------------------
            [Atoms, Bonds] = ScaleDomain(obj, Atoms, Bonds);
            
            % ---------------------------------------------------------
            % 11. Show visualization and statistics
            % ---------------------------------------------------------
            VisualizeNetwork(obj, Atoms, Bonds, Nvec);
            % VisualizeNetwork should check obj.flags.iplot internally

            % ---------------------------------------------------------
            % 12. Computes
            % ---------------------------------------------------------
            order = ComputeOrder(obj, Atoms, Bonds);

            % ---------------------------------------------------------
            % 13. Write data files
            % ---------------------------------------------------------
            WriteDataFiles(obj, Atoms, Bonds, Nvec, LDpot, order);

        end

    end


end
end