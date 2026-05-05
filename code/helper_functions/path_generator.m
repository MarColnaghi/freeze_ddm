function [path] = path_generator(varargin)
    opt = inputParser;
    addParameter(opt, 'bouts_id', 'imm2_mob2_pc4');
    addParameter(opt, 'folder', []);
    addParameter(opt, 'model', []);
    addParameter(opt, 'imfirst', false);
    parse(opt, varargin{:});

    %% 1. Dynamic Root Discovery
    % Get the full path of this script
    current_file_path = mfilename('fullpath'); 
    
    % Find the index where 'freeze_ddm' ends in the path string
    root_folder_name = 'freeze_ddm';
    idx = strfind(current_file_path, root_folder_name);
    
    if isempty(idx)
        error('Could not find the project root folder "%s" in the path.', root_folder_name);
    end
    
    % Extract everything up to 'freeze_ddm'
    project_root = current_file_path(1 : idx + length(root_folder_name) - 1);

    %% 2. Define Sub-Directories
    res_base = fullfile(project_root, 'model_results');
    fig_base = fullfile(project_root, 'figures');
    dat_base = fullfile(project_root, 'datasets');

    %% 3. Path Construction
    % Use the input parameters to build specific paths
    if ~isempty(opt.Results.folder)
        path.results = fullfile(res_base, opt.Results.folder, opt.Results.model);
        path.fig = fullfile(fig_base, opt.Results.folder, opt.Results.model);
        
        if ~isempty(opt.Results.model)
            path.filename.elbo = fullfile(path.results, 'elbo.mat');
            path.filename.est_params = fullfile(path.results, 'est_params.mat');
        end
    else
        path.results = res_base;
        path.fig = fig_base;
    end

    path.code = fullfile(project_root, 'code');
    
    % Dataset logic
    sub_folder = 'fill_in_combinations_mobfirst';
    if opt.Results.imfirst
        sub_folder = 'fill_in_combinations_imfirst';
    end

    path.exp_data = fullfile('/Users', 'marcocolnaghi', 'experimental_data', '004--social_ddm', 'dataset_1');
    path.dataset = fullfile(dat_base, sub_folder, opt.Results.bouts_id);
    path.cache_path = fullfile(dat_base, 'caches');
end