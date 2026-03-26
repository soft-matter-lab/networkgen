function [Atoms] = AddAtomsRandom(obj)
% -------------------------------------------------------------------------
% AddAtomsRandom
% - Scatter nodes in 2D with minimum spacing using a uniform grid
% - Reads all needed parameters from obj.domain
%
% OUTPUT:
%   Atoms : [ ID | X | Y | Z | ... ]
% -------------------------------------------------------------------------

    % --------- Unpack domain ---------
    xlo = obj.domain.xlo;
    xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;
    yhi = obj.domain.yhi;

    Max_atom                  = obj.domain.Max_atom;
    node_scatter_max_tries    = obj.domain.node_scatter_max_tries;
    max_tries_per_node_sample = obj.domain.max_tries_per_node_sample;
    min_node_sep2             = obj.domain.min_node_sep2;

    dmin = sqrt(min_node_sep2);

    % --------- Grid setup ---------
    Lx = xhi - xlo;
    Ly = yhi - ylo;

    h  = dmin;
    nx = max(1, ceil(Lx / h));
    ny = max(1, ceil(Ly / h));

    gridHeads = zeros(ny, nx, 'int32');
    nextIdx   = zeros(Max_atom, 1, 'int32');

    % Keep old random-network atom layout for compatibility
    Atoms  = zeros(Max_atom, 10);
    N_atom = 0;

    % --------- Helpers ---------
    function [ci,cj] = coord2cell(x, y)
        cx = floor((x - xlo) / h) + 1;
        cy = floor((y - ylo) / h) + 1;

        if cx < 1
            cx = 1;
        elseif cx > nx
            cx = nx;
        end

        if cy < 1
            cy = 1;
        elseif cy > ny
            cy = ny;
        end

        ci = cx;
        cj = cy;
    end

    function ok = passes_minsep(xi, yi)
        [ci, cj] = coord2cell(xi, yi);
        ok = true;

        for dj = -1:1
            yj = cj + dj;
            if (yj < 1) || (yj > ny)
                continue;
            end

            for di = -1:1
                xi_cell = ci + di;
                if (xi_cell < 1) || (xi_cell > nx)
                    continue;
                end

                head = gridHeads(yj, xi_cell);
                k = head;

                while k ~= 0
                    dx = xi - Atoms(k,2);
                    dy = yi - Atoms(k,3);

                    if (dx*dx + dy*dy) < min_node_sep2
                        ok = false;
                        return;
                    end

                    k = nextIdx(k);
                end
            end
        end
    end

    function insert_into_grid(idx)
        [ci, cj] = coord2cell(Atoms(idx,2), Atoms(idx,3));
        head = gridHeads(cj, ci);
        nextIdx(idx) = head;
        gridHeads(cj, ci) = int32(idx);
    end

    % --------- Main loop ---------
    global_scatter_tries = 0;
    tstart = tic;

    while (N_atom < Max_atom) && (global_scatter_tries < node_scatter_max_tries)

        global_scatter_tries = global_scatter_tries + 1;

        accepted = false;
        per_node_tries = 0;

        while (~accepted) && (per_node_tries < max_tries_per_node_sample)

            per_node_tries = per_node_tries + 1;

            xi = xlo + Lx * rand;
            yi = ylo + Ly * rand;
            zi = 0.0;

            if N_atom == 0
                accepted = true;
            else
                accepted = passes_minsep(xi, yi);
            end
        end

        if ~accepted
            continue;
        end

        N_atom = N_atom + 1;

        Atoms(N_atom,1) = N_atom;
        Atoms(N_atom,2) = xi;
        Atoms(N_atom,3) = yi;
        Atoms(N_atom,4) = zi;

        insert_into_grid(N_atom);

    end

    fprintf('   Placed %d atoms in %4.4f sec\n', N_atom, toc(tstart));

    if N_atom < Max_atom
        warning('Requested %d atoms, placed %d atoms.', Max_atom, N_atom);
        if N_atom == 0
            error('No atoms placed - aborting.');
        end
    end

    Atoms = Atoms(1:max(N_atom,1), :);

end