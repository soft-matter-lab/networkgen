function [Atoms, LatticeData] = AddAtomsHex(obj)
% -------------------------------------------------------------------------
% AddAtomsHex
% - Generate a 2D hexagonal (triangular) lattice
% - No geometric disorder applied here
%
% INPUT:
%   obj : network class object
%
% OUTPUT:
%   Atoms       : atom array
%   LatticeData : struct with idx_map, Nx, Ny
% -------------------------------------------------------------------------

    xlo = obj.domain.xlo;
    xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;
    yhi = obj.domain.yhi;

    Lx = xhi - xlo;
    Ly = yhi - ylo;

    % SetupDomain already stored the absolute min spacing
    a = obj.domain.min_node_sep;

    edgeTol = 0.25 * a;
    dy = a * sqrt(3) / 2;

    Ny_est = floor(Ly / dy) + 2;
    Nx_est = floor(Lx / a) + 3;

    idx_map = zeros(Ny_est, Nx_est);

    maxNodes = Ny_est * Nx_est;
    x_all = zeros(maxNodes,1);
    y_all = zeros(maxNodes,1);

    nat = 0;

    for iy = 1:Ny_est

        y = ylo + (iy-1) * dy;

        if (y < ylo) || (y > yhi)
            continue;
        end

        x_offset = 0;
        if mod(iy,2) == 1
            x_offset = a/2;
        end

        for ix = 1:Nx_est

            x = xlo + (ix-1)*a + x_offset;

            if (x < xlo) || (x > xhi)
                continue;
            end

            nat = nat + 1;
            x_all(nat) = x;
            y_all(nat) = y;
            idx_map(iy, ix) = nat;

        end
    end

    % Trim
    x_all = x_all(1:nat);
    y_all = y_all(1:nat);

    % Fixed boundary nodes
    isFixed = (x_all <= xlo + edgeTol) | (x_all >= xhi - edgeTol) | ...
              (y_all <= ylo + edgeTol) | (y_all >= yhi - edgeTol);

    % Keep old lattice atom layout for compatibility
    % [id x y z type isFixed]
    Atoms = zeros(nat, 10);
    Atoms(:,1) = (1:nat).';
    Atoms(:,2) = x_all;
    Atoms(:,3) = y_all;
    Atoms(:,4) = 0.0;
    Atoms(:,5) = 1;
    Atoms(:,6) = double(isFixed);

    LatticeData.idx_map = idx_map;
    LatticeData.Nx      = size(idx_map,2);
    LatticeData.Ny      = size(idx_map,1);

end