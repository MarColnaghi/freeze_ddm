function thresholds = define_thresholds

thresholds.pc = 17; thresholds.fill_in_mob = 3; thresholds.fill_in_imm = 3;
thresholds.le_window_sl = [26 45];%[20 65];
thresholds.le_window_fl = [13 32];%[5 50];

thresholds.time_window = [120 540];
thresholds.sp_window = [45 75; 100 130; 150 180; 210 240; 300 330; 390 420; -80 -50] - 30;
thresholds.freeze_dur = 1;
thresholds.mob_windows(2,:) = [80 520];