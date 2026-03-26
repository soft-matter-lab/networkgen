classdef networklog < handle
% =========================================================================
% log  -  Logging subclass for the network framework
%
% Provides two distinct logging channels:
%
%   1. CONSOLE LOG  (command-window history)
%      Any function that would normally call fprintf() should instead call
%      obj.log.print(...) with the same sprintf-style arguments.  Messages
%      are stored in an internal buffer and can later be flushed to a text
%      file that mirrors exactly what appeared in the command window.
%
%   2. NETWORK STATS LOG  (structured properties file)
%      Key network quantities (atom counts, bond counts, domain geometry,
%      density, chain statistics, order parameters, potential parameters)
%      are stored in a typed struct via record() / recordMany().
%      writeNetworkLog() formats and writes this struct to a readable file.
%
% TYPICAL USAGE IN A PIPELINE FUNCTION
% -------------------------------------------------------
%   % Instead of: fprintf('   Placed %d atoms\n', N);
%   obj.log.print('   Placed %d atoms\n', N);
%
%   % Record a statistic for the structured log
%   obj.log.record('atom_count', size(Atoms,1));
%
%   % At the end of generateNetwork, flush both files
%   obj.log.setReplicate(ii);   % stamp all output filenames with N001, N002 ...
%   obj.log.writeConsoleLog( fullfile(dir, 'console.log') );   % -> console_N001.log
%   obj.log.writeNetworkLog( fullfile(dir, 'network.log') );   % -> network_N001.log
%   obj.log.clear();          % reset messages + stats; replicate_id is kept
% -------------------------------------------------------
%
% NOTE: The class is named 'log' to match the network.m property
%   declaration  `log log`.  Inside network methods, use obj.log.print()
%   rather than the bare name to avoid shadowing MATLAB's built-in log().
% =========================================================================

properties

    % ------------------------------------------------------------------ %
    %  Console buffer
    % ------------------------------------------------------------------ %

    %%% Accumulated console messages (cell array of char strings)
    messages = {};

    %%% When true, print() also echoes to the command window immediately
    %%% Set to false to suppress all output and only buffer
    echo_to_console (1,1) logical = true;

    % ------------------------------------------------------------------ %
    %  Structured network statistics
    %
    %  stats is a struct of named fields.  Each field holds a sub-struct
    %  with fields: value, unit, fmt, description.
    %  The schema is populated by record() and written by writeNetworkLog().
    % ------------------------------------------------------------------ %

    %%% Flat key->entry map stored as a struct (dynamic field names)
    stats = struct();

    %%% Ordered list of stat keys for deterministic write order
    stat_keys = {};

    % ------------------------------------------------------------------ %
    %  Replicate tracking
    % ------------------------------------------------------------------ %

    %%% Current replicate index.  Set via setReplicate(n) at the top of
    %%% each loop iteration in generateNetwork.  Stamped into every output
    %%% filename and into the network log header.  Not reset by clear() --
    %%% it persists until explicitly changed by the next setReplicate() call.
    replicate_id = 1;
    ii = 1
    Nreplicates = 1

    seed = []
    sample_suffix = ''
    replicate_suffix = ''
    sample_label = ''

    lammps_data_file = ''
    lammps_viz_file  = ''
    bond_table_file  = ''
    log_file         = ''
    pot_file         = ''

end % properties


methods

    % =================================================================== %
    %  CONSTRUCTOR
    % =================================================================== %
    function obj = networklog()
        % Nothing to initialise beyond property defaults.
    end


    % =================================================================== %
    %  REPLICATE MANAGEMENT
    % =================================================================== %

    function setReplicate(obj, n)
    % SETREPLICATE  Set the current replicate index.
    %
    %   obj.log.setReplicate(ii)
    %
    %   Call once at the top of each iteration of the Nreplicates loop.
    %   All subsequent writeConsoleLog / writeNetworkLog / writeLogs calls
    %   will automatically insert '_N<NNN>' into the filename stem before
    %   the extension, e.g.:
    %
    %     console.log   ->  console_N001.log
    %     network.log   ->  network_N003.log
    %
    %   The replicate_id is NOT cleared by clear() so it stays correct
    %   for the duration of each replicate even after the buffers are reset.

        obj.replicate_id = max(1, round(n));
        obj.record('replicate_id', obj.replicate_id, ...
            'desc', 'Replicate index');
    end


    function clear(obj)
    % CLEAR  Reset message buffer and stats (call between replicates).
    %
    %   obj.log.clear()
    %
    %   NOTE: replicate_id is intentionally preserved so it remains valid
    %   if clear() is called mid-replicate.  Update it with setReplicate().

        obj.messages  = {};
        obj.stats     = struct();
        obj.stat_keys = {};
        % replicate_id is NOT reset here -- it persists until setReplicate()
    end


    function snapshot = snapshotStats(obj)
    % SNAPSHOTSTATS  Return a copy of the current stats struct (for archiving).
    %
    %   snap = obj.log.snapshotStats();

        snapshot = obj.stats;
    end


    % =================================================================== %
    %  CONSOLE LOGGING
    % =================================================================== %

    function print(obj, fmt, varargin)
    % PRINT  Buffer a formatted message (and optionally echo it).
    %
    %   obj.log.print('   Placed %d atoms in %.2f sec\n', N, t)
    %
    %   Follows exactly the same syntax as fprintf(fmt, ...).
    %   A trailing newline is NOT added automatically; include \n as needed.

        msg = sprintf(fmt, varargin{:});
        obj.messages{end+1} = msg;

        if obj.echo_to_console
            fprintf('%s', msg);
        end
    end


    function printSection(obj, title)
    % PRINTSECTION  Write a visible section divider to the console log.
    %
    %   obj.log.printSection('Add Bonds')
    %   ->   prints:   --- Add Bonds ---
    %
    %   Useful for visually separating pipeline stages in the log file.

        divider = sprintf('\n--- %s ---\n', title);
        obj.messages{end+1} = divider;

        if obj.echo_to_console
            fprintf('%s', divider);
        end
    end


    function writeConsoleLog(obj, filepath)
    % WRITECONSOLELOG  Write all buffered messages to a replicate-stamped file.
    %
    %   obj.log.writeConsoleLog('/path/to/console.log')
    %   -> writes to  /path/to/console_N001.log  (for replicate_id = 1)
    %
    %   The replicate suffix is inserted before the final extension so the
    %   base path can stay constant across the replicates loop.  If the
    %   path has no extension, the suffix is appended to the end.

        stamped = obj.stamp_path(filepath);

        fid = fopen(stamped, 'w');
        if fid < 0
            warning('log:writeConsoleLog:cannotOpen', ...
                    'Could not open "%s" for writing.', stamped);
            return;
        end

        for k = 1:numel(obj.messages)
            fprintf(fid, '%s', obj.messages{k});
        end

        fclose(fid);
        fprintf('   [log] Wrote console log -> %s\n', stamped);
    end


    % =================================================================== %
    %  STRUCTURED STATISTICS
    % =================================================================== %

    function record(obj, key, value, varargin)
    % RECORD  Store a single named statistic.
    %
    %   obj.log.record(key, value)
    %   obj.log.record(key, value, 'unit', unit_str)
    %   obj.log.record(key, value, 'unit', unit_str, 'desc', desc_str)
    %   obj.log.record(key, value, 'unit', unit_str, 'fmt',  fmt_str)
    %
    % PARAMETERS
    %   key   : char  - field identifier (must be a valid MATLAB identifier)
    %   value : any   - the statistic value (scalar, string, etc.)
    %
    % OPTIONAL NAME-VALUE PAIRS
    %   'unit'  : string  display units (default '')
    %   'fmt'   : string  printf format for numeric display (default '%.6g')
    %   'desc'  : string  human-readable description (default = key)
    %
    % EXAMPLES
    %   obj.log.record('atom_count',  size(Atoms,1),  'desc', 'Number of atoms')
    %   obj.log.record('rho0',        LDpot.rho0,     'unit', '1/A^2', 'fmt', '%.6f')
    %   obj.log.record('geometry',    obj.architecture.geometry)

        % Parse optional name-value pairs
        p_unit = '';
        p_fmt  = '%.6g';
        p_desc = key;

        ii = 1;
        while ii <= numel(varargin) - 1
            switch lower(varargin{ii})
                case 'unit';  p_unit = varargin{ii+1};
                case 'fmt';   p_fmt  = varargin{ii+1};
                case 'desc';  p_desc = varargin{ii+1};
            end
            ii = ii + 2;
        end

        entry.value = value;
        entry.unit  = p_unit;
        entry.fmt   = p_fmt;
        entry.desc  = p_desc;

        % Store entry and track insertion order
        if ~isfield(obj.stats, key)
            obj.stat_keys{end+1} = key;
        end
        obj.stats.(key) = entry;
    end


    function recordMany(obj, S)
    % RECORDMANY  Bulk-record a struct of plain values (no metadata).
    %
    %   obj.log.recordMany(struct('atom_count', N, 'bond_count', M, ...))
    %
    %   Each struct field is stored as a plain-value stat with no unit or
    %   description.  Useful for quickly capturing a large results struct.

        fields = fieldnames(S);
        for k = 1:numel(fields)
            f = fields{k};
            obj.record(f, S.(f));
        end
    end


    function recordNetworkStats(obj, Atoms, Bonds, Nvec, obj_network, LDpot, order)
    % RECORDNETWORKSTATS  Convenience method: auto-populate all standard
    %   network statistics from the live pipeline arrays and the network obj.
    %
    %   obj.log.recordNetworkStats(Atoms, Bonds, Nvec, obj_network, LDpot, order)
    %
    % INPUTS  (pass [] for quantities not yet computed)
    %   Atoms       : atom array  [N x (5+MaxNbr)]
    %   Bonds       : bond array  [M x 5]
    %   Nvec        : per-bond Kuhn count  [M x 1]
    %   obj_network : the network object (reads domain / architecture / flags)
    %   LDpot       : output of ConstructLDPotential (struct), or []
    %   order       : output of ComputeOrder (struct), or []

        % ---- Geometry ----
        if ~isempty(Atoms)
            obj.record('atom_count',    size(Atoms,1), ...
                'desc', 'Number of atoms');
            obj.record('atom_x_range',  [min(Atoms(:,2)), max(Atoms(:,2))], ...
                'desc', 'Atom X range [xlo xhi]');
            obj.record('atom_y_range',  [min(Atoms(:,3)), max(Atoms(:,3))], ...
                'desc', 'Atom Y range [ylo yhi]');
        end

        if ~isempty(Bonds)
            obj.record('bond_count',  size(Bonds,1), ...
                'desc', 'Number of bonds');
            obj.record('mean_L0',     mean(Bonds(:,4)), ...
                'unit', 'b-units', 'fmt', '%.4f', ...
                'desc', 'Mean equilibrium bond length');
            obj.record('std_L0',      std(Bonds(:,4)), ...
                'unit', 'b-units', 'fmt', '%.4f', ...
                'desc', 'Std dev of equilibrium bond length');
        end

        % ---- Domain ----
        if ~isempty(obj_network)
            dom = obj_network.domain;
            obj.record('domain_xlo',  dom.xlo, 'unit','b-units','fmt','%.4f', ...
                'desc','Domain xlo');
            obj.record('domain_xhi',  dom.xhi, 'unit','b-units','fmt','%.4f', ...
                'desc','Domain xhi');
            obj.record('domain_ylo',  dom.ylo, 'unit','b-units','fmt','%.4f', ...
                'desc','Domain ylo');
            obj.record('domain_yhi',  dom.yhi, 'unit','b-units','fmt','%.4f', ...
                'desc','Domain yhi');
            obj.record('domain_zlo',  dom.zlo, 'unit','b-units','fmt','%.4f', ...
                'desc','Domain zlo');
            obj.record('domain_zhi',  dom.zhi, 'unit','b-units','fmt','%.4f', ...
                'desc','Domain zhi');
            obj.record('boundary',    dom.boundary, ...
                'desc','Boundary condition');
            obj.record('geometry',    obj_network.architecture.geometry, ...
                'desc','Node placement geometry');
            obj.record('strand_mode', obj_network.architecture.strand_typology.mode, ...
                'desc','Strand topology mode');

            % Crosslink number density
            if ~isempty(Atoms)
                vol = (dom.xhi-dom.xlo) * (dom.yhi-dom.ylo) * (dom.zhi-dom.zlo);
                if vol > 0
                    obj.record('crosslink_density', size(Atoms,1)/vol, ...
                        'unit','1/b^3', 'fmt','%.6f', ...
                        'desc','Crosslink number density');
                end
            end
        end

        % ---- Chain statistics ----
        if ~isempty(Nvec) && ~isempty(Atoms)
            obj.record('mean_N_kuhn',   mean(Nvec), ...
                'fmt','%.4f', 'desc','Mean Kuhn segments per chain');
            obj.record('std_N_kuhn',    std(Nvec), ...
                'fmt','%.4f', 'desc','Std dev Kuhn segments per chain');
            obj.record('min_N_kuhn',    min(Nvec), ...
                'desc','Min Kuhn segments per chain');
            obj.record('max_N_kuhn',    max(Nvec), ...
                'desc','Max Kuhn segments per chain');
            obj.record('kuhn_per_crosslink', sum(Nvec)/max(size(Atoms,1),1), ...
                'fmt','%.4f', 'desc','Total Kuhn segments per crosslink');
        end

        % ---- Local density potential ----
        if ~isempty(LDpot)
            b = 1;
            if ~isempty(obj_network),  b = obj_network.domain.b;  end

            if isfield(LDpot,'rho0')
                obj.record('LD_rho0',    LDpot.rho0, ...
                    'fmt','%.6f', 'desc','Equilibrium local density (rho0)');
            end
            if isfield(LDpot,'R_lower')
                obj.record('LD_R_lower', LDpot.R_lower/b, ...
                    'unit','b', 'fmt','%.4f', 'desc','LD lower cutoff radius');
            end
            if isfield(LDpot,'R_upper')
                obj.record('LD_R_upper', LDpot.R_upper/b, ...
                    'unit','b', 'fmt','%.4f', 'desc','LD upper cutoff radius');
            end
            if isfield(LDpot,'rc')
                obj.record('LD_rc',      LDpot.rc/b, ...
                    'unit','b', 'fmt','%.4f', 'desc','BPM cutoff radius (rc)');
            end
        end

        % ---- Structural order ----
        if ~isempty(order)
            if isfield(order,'hex')
                if isfield(order.hex,'phi6_hexatic')
                    obj.record('order_phi6_hexatic',   order.hex.phi6_hexatic, ...
                        'fmt','%.4f', 'desc','Hexatic order parameter phi6');
                end
                if isfield(order.hex,'phi6_hexagonal')
                    obj.record('order_phi6_hexagonal', order.hex.phi6_hexagonal, ...
                        'fmt','%.4f', 'desc','Hexagonal order parameter phi6');
                end
            end
        end

    end % recordNetworkStats


    function writeNetworkLog(obj, filepath)
    % WRITENETWORKLOG  Write structured statistics to a replicate-stamped file.
    %
    %   obj.log.writeNetworkLog('/path/to/network.log')
    %   -> writes to  /path/to/network_N001.log  (for replicate_id = 1)
    %
    %   Format:
    %     <description padded to column 40> : <formatted value>  [unit]
    %
    %   Fields are written in insertion order.

        stamped = obj.stamp_path(filepath);

        fid = fopen(stamped, 'w');
        if fid < 0
            warning('log:writeNetworkLog:cannotOpen', ...
                    'Could not open "%s" for writing.', stamped);
            return;
        end

        % Header
        fprintf(fid, '==========================================================\n');
        fprintf(fid, '  NETWORK PROPERTIES LOG\n');
        fprintf(fid, '  Replicate:  %d\n',  obj.replicate_id);
        fprintf(fid, '  Generated:  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, '==========================================================\n\n');

        % Write each recorded stat
        COL_WIDTH = 38;   % description column width

        for k = 1:numel(obj.stat_keys)
            key   = obj.stat_keys{k};
            if ~isfield(obj.stats, key),  continue;  end
            entry = obj.stats.(key);

            desc_str = entry.desc;
            val      = entry.value;
            unit_str = entry.unit;
            fmt      = entry.fmt;

            val_str = format_value(val, fmt);

            pad_len = max(0, COL_WIDTH - numel(desc_str));
            pad     = repmat(' ', 1, pad_len);

            if isempty(unit_str)
                fprintf(fid, '  %s%s: %s\n', desc_str, pad, val_str);
            else
                fprintf(fid, '  %s%s: %s  [%s]\n', desc_str, pad, val_str, unit_str);
            end
        end

        fprintf(fid, '\n==========================================================\n');
        fclose(fid);
        fprintf('   [log] Wrote network log -> %s\n', stamped);
    end


    % =================================================================== %
    %  CONVENIENCE: write both logs at once
    % =================================================================== %

    function writeLogs(obj, console_path, network_path)
    % WRITELOGS  Stamp and write both log files in one call.
    %
    %   obj.log.writeLogs(console_path, network_path)
    %
    %   Both paths are stamped with the current replicate_id before writing.
    %   Typical call at the end of each generateNetwork loop iteration:
    %
    %     dir = obj.domain.write_location;
    %     obj.log.writeLogs( fullfile(dir,'console.log'), ...
    %                        fullfile(dir,'network.log') );
    %     obj.log.clear();

        obj.writeConsoleLog(console_path);
        obj.writeNetworkLog(network_path);
    end


    % =================================================================== %
    %  PRIVATE HELPERS
    % =================================================================== %

    function stamped = stamp_path(obj, filepath)
    % STAMP_PATH  Insert '_N<NNN>' before the file extension.
    %
    %   stamp_path('/out/console.log')   ->  '/out/console_N001.log'
    %   stamp_path('/out/network.log')   ->  '/out/network_N003.log'
    %   stamp_path('/out/result')        ->  '/out/result_N002'
    %
    %   The zero-padded width is 3 digits, matching up to 999 replicates.
    %   Increase the %03d format below if you need more.

        suffix = sprintf('_N%03d', obj.replicate_id);

        [folder, stem, ext] = fileparts(filepath);
        stamped = fullfile(folder, [stem, suffix, ext]);
    end


end % methods
end % classdef


% =========================================================================
%  FILE-LOCAL HELPER (not a method — lives outside the classdef block)
% =========================================================================
function s = format_value(val, fmt)
% FORMAT_VALUE  Convert a value to a display string using fmt hint.
%
%   Handles scalars, vectors, logicals, and strings/char arrays.

    if ischar(val) || isstring(val)
        s = char(val);
        return;
    end

    if islogical(val)
        if val;  s = 'true';  else;  s = 'false';  end
        return;
    end

    if isnumeric(val)
        if isscalar(val)
            s = sprintf(fmt, val);
        elseif isvector(val)
            % e.g. [xlo xhi] pair
            parts = arrayfun(@(v) sprintf(fmt, v), val, 'UniformOutput', false);
            s = ['[', strjoin(parts, '  '), ']'];
        else
            s = mat2str(val, 5);
        end
        return;
    end

    % Fallback: use MATLAB's default display
    s = evalc('disp(val)');
    s = strtrim(s);
end
