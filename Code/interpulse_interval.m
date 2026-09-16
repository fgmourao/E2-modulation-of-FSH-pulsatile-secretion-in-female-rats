function res = interpulse_interval(idx, dt, label)
%INTERPULSE_INTERVAL  Intervals between consecutive pulses.
%
%   res = INTERPULSE_INTERVAL(idx, dt, label)
%
%   idx   - cell array of pulse index vectors, one per subject (1-based)
%   dt    - sampling interval in minutes
%   label - optional name, for printing
%
%   res   - struct with all, mean, median, sd, per_subject, n, n_events
%
%   Note: intervals longer than the recording window cannot be observed, so
%   the mean is truncated at the upper end.

    if ~iscell(idx), idx = {idx}; 
    end

    if nargin < 3 || isempty(label), label = 'series'; 
    end

    ns = numel(idx);
    all_iv      = [];
    per_subject = nan(ns,1);
    n_events    = zeros(ns,1);

    for s = 1:ns
        v = sort(idx{s}(:));
        n_events(s) = numel(v);

        if numel(v) >= 2
            iv = diff(v) * dt;
            all_iv = [all_iv; iv];
            per_subject(s) = mean(iv);
        end
    end

    res = struct('all', all_iv, 'mean', mean(all_iv), 'median', median(all_iv), ...
                 'sd', std(all_iv), 'per_subject', per_subject, ...
                 'n', numel(all_iv), 'n_events', n_events);

    fprintf('--- Interpulse interval: %s ---\n', label);
    fprintf('  subjects: %d | events: %d | intervals: %d\n', ...
            ns, sum(n_events), res.n);
    if res.n == 0
        fprintf('  No subject with 2+ events.\n\n');
        return
    end
    
    fprintf('  mean %.1f min | median %.1f min | SD %.1f min | range %.0f-%.0f min\n', ...
            res.mean, res.median, res.sd, min(all_iv), max(all_iv));
    if sum(~isnan(per_subject)) > 1
        fprintf('  mean of subject means: %.1f min (SD across subjects %.1f)\n', ...
                mean(per_subject,'omitnan'), std(per_subject,'omitnan'));
    end

    fprintf('\n');

end
