function interactive(freezes, points, extra_one)
% Interactive Inverse Gaussian (Wald) PDF explorer
% Version: multiple parameter boxes + Update button

    % --- UI setup ---
    fig = uifigure('Name','Inverse Gaussian PDF Explorer','Position',[100 100 800 520]);

    ax = uiaxes(fig,'Position',[60 150 680 340]);
    hold(ax,'on')
    ax.XLabel.String = 'x';
    ax.YLabel.String = 'pdf';

    % Histogram on the UIAxes
    histogram(ax, freezes.durations_s, -1/120:1/10:11, ...
        'Normalization','probability', 'EdgeColor','none');

    % --- Parameter input fields (spread out) ---
    baseY = 90;  % vertical base
    spacingY = 35; % vertical spacing
    col1X = 150;  % first column X
    col2X = 400;  % second column X

    % Row 1
    uilabel(fig,'Text','\musm:', 'Position',[col1X baseY 45 22],'FontWeight','bold');
    betamu = uieditfield(fig,'numeric', ...
        'Limits',[-0.5 5],'Value',1.5, ...
        'Position',[col1X+55 baseY 80 22]);

    uilabel(fig,'Text','\thetafs:', 'Position',[col2X baseY 60 22],'FontWeight','bold');
    thetafs = uieditfield(fig,'numeric', ...
        'Limits',[-0.5 5],'Value',2, ...
        'Position',[col2X+65 baseY 80 22]);

    % Row 2
    uilabel(fig,'Text','\thetaln:', 'Position',[col1X baseY-spacingY 60 22],'FontWeight','bold');
    thetaln = uieditfield(fig,'numeric', ...
        'Limits',[-0.5 5],'Value',2, ...
        'Position',[col1X+65 baseY-spacingY 80 22]);

    uilabel(fig,'Text','\thetals:', 'Position',[col2X baseY-spacingY 60 22],'FontWeight','bold');
    thetals = uieditfield(fig,'numeric', ...
        'Limits',[-0.5 5],'Value',2, ...
        'Position',[col2X+65 baseY-spacingY 80 22]);

    % Row 3
    uilabel(fig,'Text','\thetab0:', 'Position',[col1X baseY-2*spacingY 60 22],'FontWeight','bold');
    thetab0 = uieditfield(fig,'numeric', ...
        'Limits',[-0.5 5],'Value',2, ...
        'Position',[col1X+65 baseY-2*spacingY 80 22]);

    uilabel(fig,'Text','ndt:', 'Position',[col2X baseY-2*spacingY 40 22],'FontWeight','bold');
    ndt = uieditfield(fig,'numeric', ...
        'Limits',[1e-2 0.5],'Value',0.2, ...
        'Position',[col2X+45 baseY-2*spacingY 80 22]);

    % Group parameters into array (update dynamically inside callback)
    params = [betamu.Value, thetafs.Value, thetaln.Value, thetals.Value, thetab0.Value, ndt.Value];

    % --- Update button ---
    uibutton(fig,'Text','Update','Position',[340 30 120 35], ...
        'ButtonPushedFcn', @(btn,evt) update(ax, ...
            [betamu.Value, thetafs.Value, thetaln.Value, thetals.Value, thetab0.Value, ndt.Value], ...
            freezes, points, extra_one));

    % --- Initial draw ---
    update(ax, params, freezes, points, extra_one);
end

function update(ax, params, freezes, points, extra_one)
    % Compute and plot model fit
    [~, f, fd] = nll_fly_ddm_newer(params, freezes, points, 'model_ed5', 'iid', 'p', extra_one);

    % Draw or update curve
    lineH = findobj(ax,'Type','Line');
    if isempty(lineH)
        plot(ax, fd, f, 'LineWidth', 1.8);
    else
        set(lineH, 'XData', fd, 'YData', f);
    end

    title(ax, 'Inverse Gaussian PDF (Updated)');
    xlim(ax, [fd(1) fd(end)]);
    ylim(ax, [0 max(f)*1.05]);
end
