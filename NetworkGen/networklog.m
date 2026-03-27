classdef networklog < handle
% =========================================================================
% networklog
%
% Logging helper for the network framework.
%
% This class provides two logging channels:
%
%   1. CONSOLE LOG
%      Buffer command-window-style messages via obj.log.print(...)
%
%   2. NETWORK STATS LOG
%      Record structured statistics via obj.log.record(...)
%      and write them to a readable text file.
%
% IMPORTANT:
%   Replicate-specific names and metadata are prepared upstream by
%   localPrepareReplicate(). This class does NOT stamp or rename files.
%   It simply writes to the exact filenames/paths it is given.
% =========================================================================

properties

    % ------------------------------------------------------------------ %
    % Console buffer
    % ------------------------------------------------------------------ %
    messages = {};

    % When true, print() also echoes to command window
    echo_to_console (1,1) logical = true;

    % ------------------------------------------------------------------ %
    % Structured stats
    % ------------------------------------------------------------------ %
    stats = struct();
    stat_keys = {};

    % ------------------------------------------------------------------ %
    % Replicate / runtime metadata
    % These are filled by localPrepareReplicate()
    % ------------------------------------------------------------------ %
    replicate_id = 1;
    ii = 1;
    Nreplicates = 1;

    seed = [];
    sample_suffix = '';
    replicate_suffix = '';
    sample_label = '';

    lammps_data_file = '';
    lammps_viz_file  = '';
    bond_table_file  = '';
    log_file         = '';
    pot_file         = '';

    console_log_file = '';
    network_log_file = '';

end


methods

    % =================================================================== %
    % Constructor
    % =================================================================== %
    function obj = networklog()
        % Property defaults are enough
    end


    % =================================================================== %
    % Replicate metadata
    % =================================================================== %
    function setReplicate(obj, n)
    % Optional lightweight helper.
    % localPrepareReplicate() is still the main source of metadata.

        obj.replicate_id = max(1, round(n));
        obj.ii = obj.replicate_id;

        obj.record('replicate_id', obj.replicate_id, ...
            'desc', 'Replicate index');
    end


    function clear(obj)
    % Clear buffered messages and recorded stats only.
    % Keep replicate metadata / filenames intact.

        obj.messages  = {};
        obj.stats     = struct();
        obj.stat_keys = {};
    end


    function snapshot = snapshotStats(obj)
        snapshot = obj.stats;
    end


    % =================================================================== %
    % Console logging
    % =================================================================== %
    function print(obj, fmt, varargin)

        msg = sprintf(fmt, varargin{:});
        obj.messages{end+1} = msg;

        if obj.echo_to_console
            fprintf('%s', msg);
        end
    end


    function printSection(obj, title)

        divider = sprintf('\n--- %s ---\n', title);
        obj.messages{end+1} = divider;

        if obj.echo_to_console
            fprintf('%s', divider);
        end
    end


    function writeConsoleLog(obj, filepath)
    % Write buffered console messages to the exact filepath provided.

        fid = fopen(filepath, 'w');
        if fid < 0
            warning('networklog:writeConsoleLog:cannotOpen', ...
                    'Could not open "%s" for writing.', filepath);
            return;
        end

        for k = 1:numel(obj.messages)
            fprintf(fid, '%s', obj.messages{k});
        end

        fclose(fid);

        if obj.echo_to_console
            fprintf('   [log] Wrote console log -> %s\n', filepath);
        end
    end


    % =================================================================== %
    % Structured statistics
    % =================================================================== %
    function record(obj, key, value, varargin)

        p_unit = '';
        p_fmt  = '%.6g';
        p_desc = key;

        ii = 1;
        while ii <= numel(varargin) - 1
            switch lower(varargin{ii})
                case 'unit'
                    p_unit = varargin{ii+1};
                case 'fmt'
                    p_fmt = varargin{ii+1};
                case 'desc'
                    p_desc = varargin{ii+1};
            end
            ii = ii + 2;
        end

        entry.value = value;
        entry.unit  = p_unit;
        entry.fmt   = p_fmt;
        entry.desc  = p_desc;

        if ~isfield(obj.stats, key)
            obj.stat_keys{end+1} = key;
        end
        obj.stats.(key) = entry;
    end


    function recordMany(obj, S)

        fields = fieldnames(S);
        for k = 1:numel(fields)
            f = fields{k};
            obj.record(f, S.(f));
        end
    end


    function recordNetworkStats(obj, Atoms, Bonds, Nvec, obj_network, LDpot, order)

        % ---- Replicate metadata ----
        obj.record('replicate_id', obj.replicate_id, ...
            'desc', 'Replicate index');

        if ~isempty(obj.seed)
            obj.record('seed', obj.seed, ...
                'desc', 'Random seed');
        end

        if ~isempty(obj.sample_label)
            obj.record('sample_label', obj.sample_label, ...
                'desc', 'Sample label');
        end

        % ---- Geometry ----
        if ~isempty(Atoms)
            obj.record('atom_count', size(Atoms,1), ...
                'desc', 'Number of atoms');

            obj.record('atom_x_range', [min(Atoms(:,2)), max(Atoms(:,2))], ...
                'desc', 'Atom X range [xlo xhi]');

            obj.record('atom_y_range', [min(Atoms(:,3)), max(Atoms(:,3))], ...
                'desc', 'Atom Y range [ylo yhi]');
        end

        if ~isempty(Bonds)
            obj.record('bond_count', size(Bonds,1), ...
                'desc', 'Number of bonds');

            obj.record('mean_L0', mean(Bonds(:,4)), ...
                'unit', 'b-units', 'fmt', '%.4f', ...
                'desc', 'Mean equilibrium bond length');

            obj.record('std_L0', std(Bonds(:,4)), ...
                'unit', 'b-units', 'fmt', '%.4f', ...
                'desc', 'Std dev of equilibrium bond length');
        end

        % ---- Domain / architecture ----
        if ~isempty(obj_network)
            dom = obj_network.domain;

            obj.record('domain_xlo', dom.xlo, ...
                'unit', 'b-units', 'fmt', '%.4f', 'desc', 'Domain xlo');
            obj.record('domain_xhi', dom.xhi, ...
                'unit', 'b-units', 'fmt', '%.4f', 'desc', 'Domain xhi');
            obj.record('domain_ylo', dom.ylo, ...
                'unit', 'b-units', 'fmt', '%.4f', 'desc', 'Domain ylo');
            obj.record('domain_yhi', dom.yhi, ...
                'unit', 'b-units', 'fmt', '%.4f', 'desc', 'Domain yhi');
            obj.record('domain_zlo', dom.zlo, ...
                'unit', 'b-units', 'fmt', '%.4f', 'desc', 'Domain zlo');
            obj.record('domain_zhi', dom.zhi, ...
                'unit', 'b-units', 'fmt', '%.4f', 'desc', 'Domain zhi');

            obj.record('boundary', dom.boundary, ...
                'desc', 'Boundary condition');

            obj.record('geometry', obj_network.architecture.geometry, ...
                'desc', 'Node placement geometry');

            obj.record('strand_mode', obj_network.architecture.strand_typology.mode, ...
                'desc', 'Strand topology mode');

            if ~isempty(Atoms)
                vol = (dom.xhi-dom.xlo) * (dom.yhi-dom.ylo) * (dom.zhi-dom.zlo);
                if vol > 0
                    obj.record('crosslink_density', size(Atoms,1)/vol, ...
                        'unit', '1/b^3', 'fmt', '%.6f', ...
                        'desc', 'Crosslink number density');
                end
            end
        end

        % ---- Chain statistics ----
        if ~isempty(Nvec) && ~isempty(Atoms)
            obj.record('mean_N_kuhn', mean(Nvec), ...
                'fmt', '%.4f', 'desc', 'Mean Kuhn segments per chain');

            obj.record('std_N_kuhn', std(Nvec), ...
                'fmt', '%.4f', 'desc', 'Std dev Kuhn segments per chain');

            obj.record('min_N_kuhn', min(Nvec), ...
                'desc', 'Min Kuhn segments per chain');

            obj.record('max_N_kuhn', max(Nvec), ...
                'desc', 'Max Kuhn segments per chain');

            obj.record('kuhn_per_crosslink', sum(Nvec)/max(size(Atoms,1),1), ...
                'fmt', '%.4f', 'desc', 'Total Kuhn segments per crosslink');
        end

        % ---- Local-density potential ----
        if ~isempty(LDpot)
            b = 1;
            if ~isempty(obj_network)
                b = obj_network.domain.b;
            end

            if isfield(LDpot, 'rho0')
                obj.record('LD_rho0', LDpot.rho0, ...
                    'fmt', '%.6f', ...
                    'desc', 'Equilibrium local density (rho0)');
            end

            if isfield(LDpot, 'R_lower')
                obj.record('LD_R_lower', LDpot.R_lower/b, ...
                    'unit', 'b', 'fmt', '%.4f', ...
                    'desc', 'LD lower cutoff radius');
            end

            if isfield(LDpot, 'R_upper')
                obj.record('LD_R_upper', LDpot.R_upper/b, ...
                    'unit', 'b', 'fmt', '%.4f', ...
                    'desc', 'LD upper cutoff radius');
            end

            if isfield(LDpot, 'rc')
                obj.record('LD_rc', LDpot.rc/b, ...
                    'unit', 'b', 'fmt', '%.4f', ...
                    'desc', 'BPM cutoff radius (rc)');
            end
        end

        % ---- Structural order ----
        if ~isempty(order)
            if isfield(order, 'hex')
                if isfield(order.hex, 'phi6_hexatic')
                    obj.record('order_phi6_hexatic', order.hex.phi6_hexatic, ...
                        'fmt', '%.4f', ...
                        'desc', 'Hexatic order parameter phi6');
                end

                if isfield(order.hex, 'phi6_hexagonal')
                    obj.record('order_phi6_hexagonal', order.hex.phi6_hexagonal, ...
                        'fmt', '%.4f', ...
                        'desc', 'Hexagonal order parameter phi6');
                end
            end
        end
    end


    function writeNetworkLog(obj, filepath)
    % Write structured stats to the exact filepath provided.

        fid = fopen(filepath, 'w');
        if fid < 0
            warning('networklog:writeNetworkLog:cannotOpen', ...
                    'Could not open "%s" for writing.', filepath);
            return;
        end

        fprintf(fid, '==========================================================\n');
        fprintf(fid, '  NETWORK PROPERTIES LOG\n');
        fprintf(fid, '  Replicate:      %d of %d\n', obj.ii, obj.Nreplicates);

        if ~isempty(obj.seed)
            fprintf(fid, '  Seed:           %d\n', obj.seed);
        end

        if ~isempty(obj.sample_suffix)
            fprintf(fid, '  Sample suffix:  %s\n', obj.sample_suffix);
        end

        if ~isempty(obj.replicate_suffix)
            fprintf(fid, '  Replicate tag:  %s\n', obj.replicate_suffix);
        end

        if ~isempty(obj.sample_label)
            fprintf(fid, '  Sample label:   %s\n', obj.sample_label);
        end

        fprintf(fid, '  Generated:      %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, '==========================================================\n\n');

        COL_WIDTH = 38;

        for k = 1:numel(obj.stat_keys)
            key = obj.stat_keys{k};

            if ~isfield(obj.stats, key)
                continue;
            end

            entry = obj.stats.(key);

            desc_str = entry.desc;
            val_str  = format_value(entry.value, entry.fmt);
            unit_str = entry.unit;

            pad_len = max(0, COL_WIDTH - numel(desc_str));
            pad = repmat(' ', 1, pad_len);

            if isempty(unit_str)
                fprintf(fid, '  %s%s: %s\n', desc_str, pad, val_str);
            else
                fprintf(fid, '  %s%s: %s  [%s]\n', desc_str, pad, val_str, unit_str);
            end
        end

        fprintf(fid, '\n==========================================================\n');
        fclose(fid);

        if obj.echo_to_console
            fprintf('   [log] Wrote network log -> %s\n', filepath);
        end
    end


    % =================================================================== %
    % Convenience: write both logs
    % =================================================================== %
    function writeLogs(obj, console_path, network_path)

        obj.writeConsoleLog(console_path);
        obj.writeNetworkLog(network_path);
    end

end
end


% =========================================================================
% File-local helper
% =========================================================================
function s = format_value(val, fmt)

    if ischar(val) || isstring(val)
        s = char(val);
        return;
    end

    if islogical(val)
        if val
            s = 'true';
        else
            s = 'false';
        end
        return;
    end

    if isnumeric(val)
        if isscalar(val)
            s = sprintf(fmt, val);
        elseif isvector(val)
            parts = arrayfun(@(v) sprintf(fmt, v), val, 'UniformOutput', false);
            s = ['[', strjoin(parts, '  '), ']'];
        else
            s = mat2str(val, 5);
        end
        return;
    end

    s = evalc('disp(val)');
    s = strtrim(s);
end