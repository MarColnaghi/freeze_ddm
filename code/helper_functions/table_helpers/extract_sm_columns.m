function [before_freeze, during_freeze] = extract_sm_columns(bouts_proc, motion_cache, varargin)
    
    opt = inputParser;
    addParameter(opt, 'chunk_dur', 30);
    addParameter(opt, 'ons', 'onsets');
    parse(opt, varargin{:});
    
    ons_vec = bouts_proc.(opt.Results.ons);
    off_vec = bouts_proc.ends - 1;
    fly_ids = bouts_proc.fly;
    chunk_dur = opt.Results.chunk_dur;
    n_trials = height(bouts_proc);
    max_dur = 630;

    % 1. Pre-calculate Cumulative Sums for each fly
    % This is much faster than slicing in a loop
    num_flies = length(motion_cache);
    cs_motion = cell(num_flies, 1);
    for f = 1:num_flies
        % Divide by 10 once here
        cs_motion{f} = cumsum(motion_cache(f)./ 10);
    end

    % 2. Vectorized calculation using the formula: sum(a:b) = CS(b) - CS(a-1)
    during_freeze = nan(n_trials, 1);
    before_freeze = nan(n_trials, 1);

    for i = 1:n_trials
        f = fly_ids(i);
        ons = ons_vec(i);
        
        % Capped "during" logic
        off_capped = min(off_vec(i), ons + max_dur - 1);
        
        % Calculate means using CS
        % During:
        during_freeze(i) = (cs_motion{f}(off_capped) - cs_motion{f}(ons-1)) / (off_capped - ons + 1);
        
        % Before:
        before_freeze(i) = (cs_motion{f}(ons - 1) - cs_motion{f}(ons - chunk_dur - 1)) / chunk_dur;
    end
end