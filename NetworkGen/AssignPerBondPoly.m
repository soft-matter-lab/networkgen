function Nvec = AssignPerBondPoly(obj, Bonds, Atoms)
% -------------------------------------------------------------------------
% AssignPerBondPoly
% - Assign per-bond Kuhn segment counts for a polydisperse network
%
% Supported modes:
%   obj.perbond.kuhn.poly.method
%       'geom'
%       'range'
%       'pmf'
%       'mono'
%
% Uses fields from:
%   obj.perbond.kuhn.poly
%
% OUTPUT:
%   Nvec : [Nbonds x 1]
% -------------------------------------------------------------------------

    natoms = size(Atoms,1);
    nbonds = size(Bonds,1);

    Nvec = zeros(nbonds,1);
    if nbonds == 0
        return;
    end

    pd = obj.perbond.kuhn.poly;
    b  = obj.domain.b;

    mode  = lower(pd.method);
    min_N = pd.min_value;

    if isempty(min_N)
        min_N = 1;
    end
    min_N = max(1, round(min_N));

    Lvec = Bonds(:,4);

    switch mode

        case 'geom'
            % N ~ (L/b)^2 with rounding policy
            raw = (Lvec ./ max(b,eps)).^2;

            switch lower(pd.rounding)
                case 'ceil'
                    Nvec = ceil(raw);
                case 'floor'
                    Nvec = floor(raw);
                otherwise
                    Nvec = round(raw);
            end

            Nvec = max(Nvec, min_N);

        case 'mono'
            % poly object can still request mono assignment
            if isfield(pd,'pmf_mean') && ~isempty(pd.pmf_mean)
                Nvec(:) = round(pd.pmf_mean);
            else
                Nvec(:) = min_N;
            end
            Nvec = max(Nvec, min_N);

        case 'range'
            % Map lengths monotonically to [target_min, target_max]
            Nlo = min(pd.target_min, pd.target_max);
            Nhi = max(pd.target_min, pd.target_max);

            if strcmpi(pd.range_method,'rank')
                [~, idx] = sort(Lvec, 'ascend');

                if nbonds == 1
                    Ntargets = (Nlo + Nhi) / 2;
                else
                    Ntargets = linspace(Nlo, Nhi, nbonds).';
                end

                Ntmp = zeros(nbonds,1);
                Ntmp(idx) = Ntargets;
                Nvec = round(Ntmp);

            else
                Lmin = min(Lvec);
                Lmax = max(Lvec);

                if Lmax == Lmin
                    Nvec = round(((Nlo + Nhi)/2) * ones(nbonds,1));
                else
                    t = (Lvec - Lmin) ./ (Lmax - Lmin);
                    Nvec = round(Nlo + t .* (Nhi - Nlo));
                end
            end

            Nvec = max(Nvec, min_N);
            Nvec = min(Nvec, max(pd.target_min, pd.target_max));

        case 'pmf'
            % Truncated geometric PMF on nu in [pmf_min, pmf_max]
            nu0   = pd.pmf_min;
            nuMax = max(pd.pmf_min, pd.pmf_max);
            K     = nuMax - nu0;

            targetMeanN = pd.pmf_mean;
            targetMeanK = max(0, targetMeanN - nu0);

            f = @(p) mean_k_of_p(p, K) - targetMeanK;

            p_lo = 1e-8;
            p_hi = 1-1e-8;

            Ps = linspace(1e-6, 1-1e-6, 200);
            Fs = zeros(size(Ps));

            for ii = 1:numel(Ps)
                Fs(ii) = f(Ps(ii));
            end

            bracket_found = false;
            a = NaN; btmp = NaN;

            for ii = 1:(numel(Ps)-1)
                if Fs(ii) == 0
                    a = Ps(ii);
                    btmp = Ps(ii);
                    bracket_found = true;
                    break;
                elseif Fs(ii) * Fs(ii+1) < 0
                    a = Ps(ii);
                    btmp = Ps(ii+1);
                    bracket_found = true;
                    break;
                end
            end

            if bracket_found
                if a == btmp
                    p_opt = a;
                else
                    p_opt = fzero(f, [a, btmp]);
                end
            else
                a_guess = max(eps, targetMeanK);
                p_opt = 1/(a_guess + 1);
                p_opt = min(max(p_opt, p_lo), p_hi);
            end

            p = min(max(p_opt, p_lo), p_hi);
            r = 1 - p;

            denom = 1 - r^(K+1);
            Pk = (p * (r .^ (0:K)).') / denom;

            exp_counts  = Pk * nbonds;
            base_counts = floor(exp_counts);
            remainder   = exp_counts - base_counts;

            assigned = sum(base_counts);
            deficit  = nbonds - assigned;

            if deficit > 0
                [~, order] = sort(remainder, 'descend');
                for t = 1:deficit
                    base_counts(order(t)) = base_counts(order(t)) + 1;
                end
            elseif deficit < 0
                [~, order] = sort(remainder, 'ascend');
                for t = 1:(-deficit)
                    jj = order(t);
                    if base_counts(jj) > 0
                        base_counts(jj) = base_counts(jj) - 1;
                    end
                end
            end

            N_list = zeros(nbonds,1);
            ptr = 1;
            for k = 0:K
                cnt = base_counts(k+1);
                if cnt <= 0
                    continue;
                end
                val = nu0 + k;
                N_list(ptr:ptr+cnt-1) = val;
                ptr = ptr + cnt;
            end

            if ptr <= nbonds
                N_list(ptr:nbonds) = nu0 + K;
            end

            align = 'ascend';
            if isfield(pd,'align_to_length') && ~isempty(pd.align_to_length)
                align = lower(pd.align_to_length);
            end

            if strcmp(align,'ascend')
                [~, idxL] = sort(Lvec, 'ascend');
                Nvec(idxL) = N_list;
            else
                rp = randperm(nbonds).';
                Nvec(rp) = N_list;
            end

            Nvec = max(Nvec, min_N);

        otherwise
            error('AssignPerBondPoly: unknown poly.method "%s".', pd.method);
    end

    obj.log.print('   Kuhn assignment mode: poly (%s)\n', mode);
    obj.log.print('   Kuhn-to-crosslinker ratio %0.4f\n', sum(Nvec)/max(natoms,1));
    obj.log.print('   Average chain length %0.4f\n', sum(Nvec)/max(length(Nvec),1));

end


% ===== helper =====
function mk = mean_k_of_p(p, K)
% Truncated geometric on k = 0..K with Pk ~ p (1-p)^k.
% Returns mean(k).

    r = 1 - p;
    num = (1-p) .* (1 - (K+1)*r.^K + K*r.^(K+1)) ./ p;
    den = 1 - r.^(K+1);
    mk = num ./ den;

end