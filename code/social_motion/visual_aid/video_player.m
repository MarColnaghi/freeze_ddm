% 1. Load Video and Data
[videoFile, videoPath] = uigetfile({'*.mp4;*.mov;*.avi', 'Video Files'});
videoObj = VideoReader(fullfile(videoPath, videoFile));
coords = readmatrix('MH-1CS1NorpA3LC6Chrimson_5F-20LoomOpto100mW5Hz2min33cm-FH1_1-2021-04-09T11_38_31.trajectories.csv'); 
tsData = readtable('1CS1NorpA3LC6ChR_5F-25cm-FH1_2-2021-04-09T11_38_31.fly1_scored_frames.csv');

% --- Settings ---
windowSize = 300; % How many frames to show at once in the sliding window

% 2. Setup UI Components
fig = uifigure('Name', 'Sliding Sync Player', 'Position', [50 100 1200 600]);

% Left Side: Video
axVideo = uiaxes(fig, 'Position', [50 150 550 400]);

% Right Side: Time Series (Stacked)
axMotion = uiaxes(fig, 'Position', [650 360 500 180]);
axFreezing = uiaxes(fig, 'Position', [650 150 500 180]);

% Plot the full data once
plot(axMotion, tsData.sum_motion, 'Color', [0 0.447 0.741], 'LineWidth', 1);
ylabel(axMotion, 'Sum Motion');
hold(axMotion, 'on');

plot(axFreezing, tsData.freeze_bout, 'Color', [0.85 0.325 0.098], 'LineWidth', 1);
ylabel(axFreezing, 'Freezing Bout');
xlabel(axFreezing, 'Frame');
hold(axFreezing, 'on');

% Add Center Reference Lines (Fixed at the "Current Frame" indicator)
% These lines don't move; the graph underneath them does.
lineM = xline(axMotion, 1, 'Color', 'r', 'LineWidth', 1.5, 'Label', 'Current');
lineF = xline(axFreezing, 1, 'Color', 'r', 'LineWidth', 1.5);

% Slider
sld = uislider(fig, 'Position', [100 80 1000 3], ...
    'Limits', [1, videoObj.NumFrames], ...
    'ValueChangedFcn', @(sld, event) updateFrame(sld, videoObj, coords, axVideo, axMotion, axFreezing, lineM, lineF, windowSize));

% Play Button
btn = uibutton(fig, 'state', 'Text', 'Play', 'Position', [550 30 100 30]);

% 3. Timer Setup 
% Determine a safe FrameRate
safeFPS = videoObj.FrameRate;
if isempty(safeFPS) || safeFPS <= 0
    safeFPS = 30; % Default to 30 FPS if metadata is missing
end

% Calculate Period (must be at least 0.001)
timerPeriod = max(0.001, round(1/safeFPS));

t = timer('ExecutionMode', 'fixedRate', ...
          'Period', timerPeriod, ...
          'TimerFcn', @(~,~) stepForward(sld, videoObj, coords, axVideo, axMotion, axFreezing, lineM, lineF, windowSize, btn));

btn.ValueChangedFcn = @(btn, event) togglePlayback(btn, t);
fig.CloseRequestFcn = @(src, event) deleteTimer(t, src);

% Initial Display
updateFrame(sld, videoObj, coords, axVideo, axMotion, axFreezing, lineM, lineF, windowSize);

% --- Helper Functions ---

function togglePlayback(btn, t)
    if btn.Value, btn.Text = 'Stop'; start(t);
    else, btn.Text = 'Play'; stop(t); end
end

function stepForward(sld, videoObj, coords, axV, axM, axF, lineM, lineF, win, btn)
    if sld.Value < sld.Limits(2)
        sld.Value = sld.Value + 1;
        updateFrame(sld, videoObj, coords, axV, axM, axF, lineM, lineF, win);
    else
        btn.Value = 0;
        togglePlayback(btn, timerfind('Running', 'on'));
    end
end

function updateFrame(sld, videoObj, coords, axV, axM, axF, lineM, lineF, win)
    frameIdx = round(sld.Value);
    
    % 1. Update Video
    videoObj.CurrentTime = (frameIdx - 1) / videoObj.FrameRate;
    if hasFrame(videoObj)
        frame = readFrame(videoObj);
        frameCoords = coords(frameIdx, :);
        points = reshape(frameCoords, [2, 5])'; 
        radii = repmat(15, 5, 1);
        annotatedFrame = insertShape(frame, 'Circle', [points, radii], 'LineWidth', 3, 'Color', 'cyan');
        imshow(annotatedFrame, 'Parent', axV);
        title(axV, sprintf('Frame: %d', frameIdx));
    end
    
    % 2. Update Sliding Window (The "Magic" part)
    % We center the X-axis on frameIdx
    halfWin = win / 2;
    newXLim = [frameIdx - halfWin, frameIdx + halfWin];
    
    axM.XLim = newXLim;
    axF.XLim = newXLim;
    
    % Move the red line to follow the frame index
    lineM.Value = frameIdx;
    lineF.Value = frameIdx;
end

function deleteTimer(t, fig)
    stop(t); delete(t); delete(fig);
end