col = cmapper();
thresholds = define_thresholds();
imfirst = false;

for threshold_imm = [3]
   % for threshold_mob =  11:-2:1
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

                [bouts] = load_flies_new(thresholds, 'paths', paths, 'save', 'bouts', 'edit_filename', false, 'imfirst', imfirst);
            end
        end
  %  end
end