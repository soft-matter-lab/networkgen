function LDpot = ConstructLDPotential(obj, Atoms, Bonds, Nvec)
% -------------------------------------------------------------------------
% ConstructLDPotential
% - Construct local-density potential parameters and table
% - Reads settings from obj.domain and obj.pot
%
% INPUT:
%   obj   : network object
%   Atoms : atom array
%   Bonds : bond array
%   Nvec  : per-bond Kuhn segment counts
%
% OUTPUT:
%   LDpot : struct with local-density potential parameters
% -------------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % Basic checks
    % ---------------------------------------------------------------------
    if nargin < 4
        error('ConstructLDPotential: requires obj, Atoms, Bonds, and Nvec.');
    end

    Atom_count = size(Atoms,1);
    Bond_count = size(Bonds,1); %#ok<NASGU>

    if Atom_count == 0
        error('ConstructLDPotential: Atoms is empty.');
    end

    if isempty(Nvec)
        error('ConstructLDPotential: Nvec is empty.');
    end

    % ---------------------------------------------------------------------
    % Unpack domain / potential settings
    % ---------------------------------------------------------------------
    xlo = obj.domain.xlo;
    xhi = obj.domain.xhi;
    ylo = obj.domain.ylo;
    yhi = obj.domain.yhi;

    b = obj.domain.b;

    kLD     = obj.pot.k_LD;
    N_rho   = obj.pot.N_rho;
    rho_min = obj.pot.rho_min;
    rho_max = obj.pot.rho_max;

    if isempty(kLD)
        kLD = 0.414;
    end
    if isempty(N_rho)
        N_rho = 100000;
    end
    if isempty(rho_min)
        rho_min = 0.0;
    end
    if isempty(rho_max)
        rho_max = 500.0;
    end

    if N_rho < 2
        error('ConstructLDPotential: obj.pot.N_rho must be >= 2.');
    end

    drho = (rho_max - rho_min) / (N_rho - 1);

    % ---------------------------------------------------------------------
    % Derived network quantities
    % ---------------------------------------------------------------------
    Total_kuhn_segment = sum(Nvec);

    if Total_kuhn_segment <= 0
        error('ConstructLDPotential: sum(Nvec) must be positive.');
    end

    sig_c = 0.5 * b * sqrt(Total_kuhn_segment / Atom_count);

    atom_area = (xhi - xlo) * (yhi - ylo);
    if atom_area <= 0
        error('ConstructLDPotential: non-positive in-plane domain area.');
    end

    atom_density = Atom_count / atom_area; %#ok<NASGU>

    % ---------------------------------------------------------------------
    % Construct local-density potential parameters
    % ---------------------------------------------------------------------
    R2 = 4.0 * sig_c;
    rho0 = 0.8 * (R2 / sig_c)^2;

    R1 = 0.8 * sig_c;
    rc = 2.0 * sig_c;

    % Density vector and harmonic density potential
    rho_vec = linspace(rho_min, rho_max, N_rho + 1).';
    pot_density = kLD * (rho_vec - rho0).^2;

    % ---------------------------------------------------------------------
    % Pack output struct
    % ---------------------------------------------------------------------
    LDpot = struct();

    LDpot.N_LD        = 1;
    LDpot.N_rho       = N_rho;
    LDpot.R_lower     = R1;
    LDpot.R_upper     = R2;
    LDpot.rc          = rc;
    LDpot.rho_min     = rho_min;
    LDpot.rho0        = rho0;
    LDpot.rho_max     = rho_max;
    LDpot.drho        = drho;
    LDpot.pot_density = pot_density;

    % Optional extra bookkeeping
    LDpot.sig_c = sig_c;

    % ---------------------------------------------------------------------
    % Logging
    % ---------------------------------------------------------------------
    obj.log.print('   Constructed LD potential with parameters:\n');
    obj.log.print('   Target equilibrium density rho0 = %.4f\n', rho0);
    obj.log.print('   Lower cutoff R1 = %.4f * b\n', R1 / b);
    obj.log.print('   Upper cutoff R2 = %.4f * b\n', R2 / b);
    obj.log.print('   BPM/spring repulsion cutoff rc = %.4f * b\n', rc / b);
end