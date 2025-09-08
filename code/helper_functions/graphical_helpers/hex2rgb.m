function rgb = hex2rgb(hex)
    % Converts a hexadecimal color string to an RGB triplet.
    %
    % Input:
    %   hex - A string representing the hexadecimal color (e.g., '#FF5733' or 'FF5733').
    %
    % Output:
    %   rgb - A 1x3 array representing the RGB triplet (values between 0 and 1).

    % Remove '#' if it exists
    if startsWith(hex, '#')
        hex = hex(2:end);
    end

    % Validate the input
    if length(hex) ~= 6 || ~all(ismember(hex, '0123456789ABCDEFabcdef'))
        error('Input must be a valid 6-character hexadecimal string.');
    end

    % Convert hexadecimal to RGB
    r = hex2dec(hex(1:2)) / 255;
    g = hex2dec(hex(3:4)) / 255;
    b = hex2dec(hex(5:6)) / 255;

    % Combine into an RGB triplet
    rgb = [r, g, b];
end
