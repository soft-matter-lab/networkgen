function [] = localPrepareReplicate(obj, ii, Nreplicates)
% -------------------------------------------------------------------------
% localPrepareReplicate
% - Prepare replicate-specific runtime info
% - Sets RNG seed
% - Builds sample / replicate suffix strings
% - Builds replicate-specific output filenames
%
% INPUT:
%   obj         : network object
%   ii          : current replicate index
%   Nreplicates : total number of replicates
%
% OUTPUT:
%   none
%   (writes replicate-specific runtime info into obj.log)
% -------------------------------------------------------------------------

    obj.log.print('Generating network replicate %d of %d...\n', ii, Nreplicates);

    %% ---------------------- 1. Choose replicate seed --------------------
    if obj.flags.imanualseed
        seed_in = obj.domain.seed;

        if isempty(seed_in)
            error('localPrepareReplicate: manual seed mode is on, but obj.domain.seed is empty.');
        end

        if numel(seed_in) == 1
            seed_now = seed_in;
        else
            if numel(seed_in) < Nreplicates
                error(['localPrepareReplicate: not enough manual seeds for the ' ...
                       'number of replicates. Expected %d, got %d.'], ...
                       Nreplicates, numel(seed_in));
            end
            seed_now = seed_in(ii);
        end
    else
        seed_now = randi(1e6);
    end

    rng(seed_now);
    obj.log.seed = seed_now;

    %% ---------------------- 2. Sample suffix ----------------------------
    smp_number = obj.domain.smp_number;
    sample_suffix = sprintf('SMP%04d', smp_number);

    %% ---------------------- 3. Replicate suffix -------------------------
    if Nreplicates > 1
        replicate_suffix = sprintf('N%04d', ii);
    else
        replicate_suffix = sprintf('N%04d', 1);
    end

    %% ---------------------- 4. Network-type prefix ----------------------
    network_prefix = resolveNetworkTypePrefix(obj);

    %% ---------------------- 5. Build label ------------------------------
    % Example: PD_SMP0001_N0003
    sample_label = sprintf('%s_%s_%s', ...
        network_prefix, sample_suffix, replicate_suffix);

    %% ---------------------- 6. Store runtime info -----------------------
    obj.log.ii               = ii;
    obj.log.Nreplicates      = Nreplicates;
    obj.log.sample_suffix    = sample_suffix;
    obj.log.replicate_suffix = replicate_suffix;
    obj.log.sample_label     = sample_label;

    %% ---------------------- 7. Resolve output filenames -----------------
    % Keep obj.domain.* as base names, write resolved replicate-specific
    % names into obj.log so we do not permanently mutate the base names.

    if obj.flags.savemode
        obj.log.lammps_data_file = sprintf('%s_%s.dat', ...
            obj.domain.lammps_data_file, sample_label);

        obj.log.lammps_viz_file = sprintf('%s_%s.dat', ...
            obj.domain.lammps_viz_file, sample_label);

        obj.log.bond_table_file = sprintf('%s_%s.table', ...
            obj.domain.bond_table_file, sample_label);

        obj.log.log_file = sprintf('%s.log', sample_label);

        obj.log.pot_file = sprintf('%s.localdensity.table', sample_label);

        obj.log.console_log_file = sprintf('console_%s.log', sample_label);
        obj.log.network_log_file = sprintf('network_%s.log', sample_label);

    else
        % Manual / raw naming mode: keep base names
        obj.log.lammps_data_file = sprintf('%s.dat', obj.domain.lammps_data_file);
        obj.log.lammps_viz_file  = sprintf('%s.dat', obj.domain.lammps_viz_file);
        obj.log.bond_table_file  = sprintf('%s.table', obj.domain.bond_table_file);
        obj.log.log_file         = sprintf('%s.log', obj.domain.lammps_data_file);
        obj.log.pot_file         = sprintf('%s.localdensity.table', obj.domain.bond_table_file);
        obj.log.console_log_file = 'console.log';
        obj.log.network_log_file = 'network.log';
    end

end