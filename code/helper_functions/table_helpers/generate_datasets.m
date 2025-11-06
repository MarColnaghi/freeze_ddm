clear all
col = cmapper();
thresholds = define_thresholds();
imfirst = true;

for threshold_imm = 4:-1:3
    % for threshold_mob =  4:-1:1
        for threshold_pc = 4

            close all

            id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm - 1, threshold_imm - 1, threshold_pc);
            paths = path_generator('bouts_id', id_code, 'imfirst', imfirst);
            mkdir(paths.dataset)

            thresholds.fill_in_imm = threshold_imm;
            thresholds.fill_in_mob = threshold_imm;
            thresholds.pc = threshold_pc;

            if isfile(fullfile(paths.dataset, 'bouts.mat'))
                disp('already have this dataset')
            else
                if imfirst
                    [bouts] = load_flies_new(thresholds, 'paths', paths, 'save', 'bouts', 'edit_filename', false, 'imfirst', imfirst);
                else
                    [bouts] = load_flies_new(thresholds, 'paths', paths, 'save', 'bouts', 'edit_filename', false, 'imfirst', imfirst);
                end
            end
        end
        %end
end