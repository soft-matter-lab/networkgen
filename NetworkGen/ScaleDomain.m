function [Atoms, Bonds] = ScaleDomain(obj, Atoms, Bonds)
% -------------------------------------------------------------------------
% ScaleDomain
% - Uniformly rescale domain bounds, atom coordinates, and bond lengths
% - Reads scale factor from obj.domain.scale
%
% INPUT:
%   obj   : network object
%   Atoms : atom array [ID x y z ...]
%   Bonds : bond array [ID i j L type]
%
% OUTPUT:
%   Atoms : scaled atom array
%   Bonds : scaled bond array
%
% SIDE EFFECT:
%   Updates obj.domain.xlo/xhi/ylo/yhi/zlo/zhi in place
% -------------------------------------------------------------------------

    if nargin < 3
        error('ScaleDomain: requires obj, Atoms, and Bonds.');
    end

    scaleFactor = obj.domain.scale;

    if isempty(scaleFactor)
        scaleFactor = 1.0;
    end

    if ~isscalar(scaleFactor) || ~isnumeric(scaleFactor) || ~isfinite(scaleFactor)
        error('ScaleDomain: obj.domain.scale must be a finite numeric scalar.');
    end

    if scaleFactor <= 0
        error('ScaleDomain: obj.domain.scale must be positive.');
    end

    % No-op case
    if scaleFactor == 1
        obj.log.print('   Domain scaling skipped (scale = 1).\n');
        return;
    end

    obj.log.print('   Scaling domain by factor %.6g\n', scaleFactor);

    % ---------------------------------------------------------------------
    % Scale domain bounds
    % ---------------------------------------------------------------------
    obj.domain.xlo = obj.domain.xlo * scaleFactor;
    obj.domain.xhi = obj.domain.xhi * scaleFactor;

    obj.domain.ylo = obj.domain.ylo * scaleFactor;
    obj.domain.yhi = obj.domain.yhi * scaleFactor;

    obj.domain.zlo = obj.domain.zlo * scaleFactor;
    obj.domain.zhi = obj.domain.zhi * scaleFactor;

    % ---------------------------------------------------------------------
    % Scale atom coordinates
    % ---------------------------------------------------------------------
    if ~isempty(Atoms) && size(Atoms,2) >= 4
        Atoms(:,2:4) = Atoms(:,2:4) * scaleFactor;
    end

    % ---------------------------------------------------------------------
    % Scale bond equilibrium lengths
    % ---------------------------------------------------------------------
    if ~isempty(Bonds) && size(Bonds,2) >= 4
        Bonds(:,4) = Bonds(:,4) * scaleFactor;
    end

    obj.log.print('   Scaled domain bounds and network geometry.\n');
    obj.log.print('   New x range: [%.6g, %.6g]\n', obj.domain.xlo, obj.domain.xhi);
    obj.log.print('   New y range: [%.6g, %.6g]\n', obj.domain.ylo, obj.domain.yhi);
    obj.log.print('   New z range: [%.6g, %.6g]\n', obj.domain.zlo, obj.domain.zhi);

end