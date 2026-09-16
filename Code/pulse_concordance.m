function res = pulse_concordance(idxA, idxB, N, window, lags, nsim, dt, ids, seed)

%pulse_concordance  Mean SC and mean SC max across subjects.
%
%   res = pulse_concordance(idxA, idxB, N, window, lags, nsim, dt)
%   res = pulse_concordance(idxA, idxB, N, window, lags, nsim, dt, ids, seed)
%
%   Specific Concordance index of Guardabasso et al., Acta Endocrinol
%   1991;124:208-218. Takes pulse indices directly, averages the SC spectrum
%   across subjects, and tests the maximum against a Monte Carlo null that
%   reproduces the design: same subjects, same event counts, same N.
%
%   idxA, idxB - cell arrays of pulse index vectors (1-based), one per
%                subject, in the same subject order
%   N          - samples per series
%   window     - time window in observations (odd; 3 in the original paper)
%   lags       - lags tested, in samples (e.g. -8:8)
%   nsim       - Monte Carlo simulations
%   dt         - sampling interval in minutes
%   ids        - optional subject labels
%   seed       - optional RNG seed; without it percentiles vary between runs
%
%   res        - struct with spectra, mean_spectrum, mean_SCmax, lag_min,
%                p, p95, p99, null_baseline, n_subjects, events
%
%   Conventions: positive lag means an event in A preceded one in B. When
%   sliding, N is held constant and events of B falling outside 1..N are
%   discarded, with nB recounted among survivors. Each event is paired at
%   most once.
%
%   MEAN SC VALUES ARE NOT COMPARABLE ACROSS GROUPS. The expected SC under
%   the null is not zero: it becomes more negative as more events are
%   detected, since chance concordance grows with nA*nB/N^2 while
%   coincidences grow linearly. Compare significance across groups, not SC
%   magnitude. res.null_baseline reports the effective zero for this group.
%
%   The p95/p99 cutoffs are the null distribution of the MAXIMUM over the
%   lag range. They already correct for multiple comparisons and answer one
%   question: is there non-random concordance at any lag? Do not read the
%   spectrum point by point against them, and do not narrow the lag range
%   after seeing the result

    ns = numel(idxA);

    if numel(idxB) ~= ns
        error('idxA has %d subjects, idxB has %d.', ns, numel(idxB));
    end

    if nargin >= 9 && ~isempty(seed), rng(seed); 
    end

    if nargin < 8 || isempty(ids)
        ids = arrayfun(@(i) sprintf('subject %d', i), 1:ns, 'UniformOutput', false);
    end

    if mod(window,2) == 0
        error('window must be odd; got %d.', window);
    end

    if nsim < 10000
        warning('pulse_concordance:nsim', ...
            ['nsim = %d. The 99th percentile rests on ~%d simulations and ' ...
             'varies between runs. Use nsim >= 20000 for final results.'], ...
             nsim, max(1, round(nsim*0.01)));
    end

    % validate and store indices
    A = cell(1,ns);  B = cell(1,ns);
    for s = 1:ns
        A{s} = local_check(idxA{s}, N, ids{s}, 'A');
        B{s} = local_check(idxB{s}, N, ids{s}, 'B');
    end

    nev = zeros(ns,2);
    spectra = nan(ns, numel(lags));
    for s = 1:ns
        nev(s,:) = [numel(A{s}) numel(B{s})];
        
        for k = 1:numel(lags)
            spectra(s,k) = local_sc(A{s}, B{s}, N, lags(k), window);
        end

    end

    mean_spectrum = mean(spectra, 1, 'omitnan');
    [mean_SCmax, imax] = max(mean_spectrum);

    % Monte Carlo null

    nullmax = zeros(nsim,1);
    nullspec = zeros(nsim, numel(lags));

    for t = 1:nsim
        sim = nan(ns, numel(lags));

        for s = 1:ns
            sa = local_random(N, nev(s,1));
            sb = local_random(N, nev(s,2));

            for k = 1:numel(lags)
                sim(s,k) = local_sc(sa, sb, N, lags(k), window);
            end

        end

        nullspec(t,:) = mean(sim, 1, 'omitnan');
        nullmax(t)    = max(nullspec(t,:));

    end

    res = struct();
    res.ids           = ids;
    res.spectra       = spectra;
    res.mean_spectrum = mean_spectrum;
    res.lags          = lags;
    res.dt            = dt;
    res.mean_SCmax    = mean_SCmax;
    res.lag_min       = lags(imax) * dt;
    res.p             = mean(nullmax >= mean_SCmax);
    res.p95           = local_prctile(nullmax, 95);
    res.p99           = local_prctile(nullmax, 99);
    res.null_baseline = mean(nullspec(:));
    res.n_subjects    = ns;
    res.events        = nev;

    fprintf('--- Mean SC (%d subjects) ---\n', ns);
    fprintf('  subject        nA   nB    SCmax    lag(min)\n');

    for s = 1:ns
        [sm, im] = max(spectra(s,:));
        fprintf('  %-12s  %3d  %3d   %+.4f      %+4d\n', ...
                ids{s}, nev(s,1), nev(s,2), sm, lags(im)*dt);
    end

    fprintf('\n  Mean SCmax = %+.4f at lag %+d min\n', res.mean_SCmax, res.lag_min);
    fprintf('  H0: p95 = %.4f | p99 = %.4f  ->  p = %.4f\n', ...
            res.p95, res.p99, res.p);
    fprintf('  null baseline = %+.4f  (total %d events in A, %d in B)\n', ...
            res.null_baseline, sum(nev(:,1)), sum(nev(:,2)));

    if res.mean_SCmax > res.p99
        fprintf('  NON-RANDOM CONCORDANCE (p < 0.01)\n\n');

    elseif res.mean_SCmax > res.p95
        fprintf('  NON-RANDOM CONCORDANCE (p < 0.05)\n\n');

    else
        fprintf('  Not significant.\n\n');
    end

end

% ======================================================================
function v = local_check(v, N, id, which)

% Validate a pulse index vector

    v = v(:)';

    if isempty(v), return; end

    if any(mod(v,1) ~= 0)
        error('%s, series %s: non-integer index. If the detector reported time in minutes, convert with idx = time/dt + 1.', id, which);
    end

    if any(v < 1)
        error('%s, series %s: index < 1. The detector probably counts from 0; add 1.', id, which);
    end

    if any(v > N)
        error('%s, series %s: index %d exceeds N = %d.', id, which, max(v), N);
    end

    if numel(unique(v)) < numel(v)
        warning('%s, series %s: duplicate indices, counted once.', id, which);
        v = unique(v);
    end

    if any(diff(sort(v)) == 1)
        warning(['%s, series %s: events in adjacent samples. The null model ' ...
                 'forbids this in simulated series.'], id, which);
    end

    v = sort(v);

end

% ======================================================================
function sc = local_sc(evA, evB, N, lag, w)

% Specific Concordance for one subject at one lag

    if abs(lag) >= N, sc = NaN; return; 
    end

    h  = floor((w - 1) / 2);
    Bs = evB - lag;
    Bs = Bs(Bs >= 1 & Bs <= N);

    used = false(1, numel(Bs));
    a = 0;
    for x = evA
        for k = 1:numel(Bs)
            if ~used(k) && abs(Bs(k) - x) <= h
                used(k) = true;  a = a + 1;  break
            end
        end
    end

    fc  = (numel(evA)/N) * (numel(Bs)/N);
    fcc = 1 - (1 - fc)^w;        % generalises eq. [2]; w = 3 gives 3fc-3fc^2+fc^3
    sc  = a/N - fcc;
end

% ======================================================================
function v = local_random(N, nev)

% Random event positions, no two adjacent

    if nev == 0, v = []; return; 
    end

    for k = 1:500
        p = sort(randperm(N, nev));
        if nev < 2 || all(diff(p) > 1), v = p; return; end
    end

    v = 1:2:min(N, 2*nev);
    v = v(1:min(nev, numel(v)));

end

% ======================================================================
function q = local_prctile(x, pc)
    x = sort(x(:));  n = numel(x);
    i = (pc/100)*(n-1) + 1;
    lo = floor(i);  hi = ceil(i);  f = i - lo;
    q = x(lo)*(1-f) + x(hi)*f;
end
