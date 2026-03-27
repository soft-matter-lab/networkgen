function order = ComputeOrder(obj, Atoms, Bonds)
% -------------------------------------------------------------------------
% ComputeOrder
% - Compute structural order parameters for the generated network
% - Reads Atoms and Bonds and returns a packed order struct
%
% INPUT:
%   obj   : network object
%   Atoms : atom array
%   Bonds : bond array
%
% OUTPUT:
%   order : struct containing order-parameter results
% -------------------------------------------------------------------------

    %#ok<INUSD>
    order = struct();

    % ---------------------------------------------------------------------
    % Hexagonal order parameter
    % ---------------------------------------------------------------------
    [phi6k, phi6_hexatic, phi6_hexagonal] = ComputeHexOrder(Atoms, Bonds);

    order.hex.phi6k           = phi6k;
    order.hex.phi6_hexatic    = phi6_hexatic;
    order.hex.phi6_hexagonal  = phi6_hexagonal;

    % ---------------------------------------------------------------------
    % Logging
    % ---------------------------------------------------------------------
    obj.log.print('   Computed structural order parameters:\n');
    obj.log.print('   Hexatic order phi6 = %.4f\n', phi6_hexatic);
    obj.log.print('   Hexagonal order phi6 = %.4f\n', phi6_hexagonal);

end


function [phi6k, phi6_hexatic, phi6_hexagonal] = ComputeHexOrder(Atoms, Bonds)
% -------------------------------------------------------------------------
% ComputeHexOrder
% - Compute local and global hexagonal/hexatic order parameters from
%   network connectivity
%
% INPUT:
%   Atoms : atom array
%   Bonds : bond array [bondID id1 id2 ...]
%
% OUTPUT:
%   phi6k          : complex local hexatic order parameter per atom
%   phi6_hexatic   : global hexatic order parameter
%   phi6_hexagonal : average local hexagonal order magnitude
% -------------------------------------------------------------------------

    N = size(Atoms,1);

    phi6k    = zeros(N,1);
    phi6norm = zeros(N,1);

    if N == 0
        phi6_hexatic = 0;
        phi6_hexagonal = 0;
        return;
    end

    % ---------------------------------------------------------------------
    % Build neighbor list from bond connectivity
    % ---------------------------------------------------------------------
    neighborList = cell(N,1);

    for iB = 1:size(Bonds,1)
        a1 = Bonds(iB,2);
        a2 = Bonds(iB,3);

        if a1 >= 1 && a1 <= N && a2 >= 1 && a2 <= N
            neighborList{a1} = [neighborList{a1}, a2]; %#ok<AGROW>
            neighborList{a2} = [neighborList{a2}, a1]; %#ok<AGROW>
        end
    end

    % ---------------------------------------------------------------------
    % Compute local hexatic order parameter for each atom
    % ---------------------------------------------------------------------
    for k = 1:N

        neighbors = neighborList{k};
        numNeighbors = numel(neighbors);

        if numNeighbors == 0
            phi6k(k) = 0;
            phi6norm(k) = 0;
            continue;
        end

        angles = zeros(numNeighbors,1);

        for j = 1:numNeighbors
            vec = Atoms(neighbors(j), 2:3) - Atoms(k, 2:3);
            angles(j) = atan2(vec(2), vec(1));
        end

        phi6k(k) = sum(exp(1i * 6 * angles)) / numNeighbors;
        phi6norm(k) = abs(phi6k(k));
    end

    % ---------------------------------------------------------------------
    % Global order parameters
    % ---------------------------------------------------------------------
    phi6_hexatic   = abs(sum(phi6k) / N);
    phi6_hexagonal = sum(phi6norm) / N;

end