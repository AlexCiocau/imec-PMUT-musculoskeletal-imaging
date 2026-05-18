% % 1. Load the Data
% filename = 'C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\YZ_focused_detailed\Scan_2D_UMSmap.txt';
% % The file has 4 header lines before the actual number matrix starts
% rawData = readmatrix(filename, 'NumHeaderLines', 4);
% 
% % 2. Extract Axes and Pressure Matrix
% % Y-axis coordinates are in the first row, columns 2 through 13
% Y_raw = rawData(1, 2:end); 
% 
% % Z-axis (depth) coordinates are in the first column, rows 2 through 31
% Z_raw = rawData(2:end, 1); 
% 
% % The pressure data is the remaining 30x12 block
% P_raw = rawData(2:end, 2:end); 
% 
% % 3. Coordinate Adjustments
% % Center the Y-axis around 0 mm to clearly see the beam spread
% Y = Y_raw - mean(Y_raw);
% 
% % Shift the Z-axis so the scan starts exactly at 3 mm depth
% Z = Z_raw - Z_raw(1) + 3;
% 
% % 4. Interpolation for Smoothing
% % Create grids for the original 30x12 data resolution
% [Y_grid, Z_grid] = meshgrid(Y, Z);
% 
% % Create a much finer grid for a smooth heatmap (e.g., 200x200 resolution)
% Y_fine = linspace(min(Y), max(Y), 200);
% Z_fine = linspace(min(Z), max(Z), 200);
% [Y_grid_fine, Z_grid_fine] = meshgrid(Y_fine, Z_fine);
% 
% % Perform 2D spline interpolation to smooth the jagged edges
% P_fine = interp2(Y_grid, Z_grid, P_raw, Y_grid_fine, Z_grid_fine, 'spline');
% 
% % 5. Plotting the Heatmap
% figure('Color', 'w', 'Name', 'PMUT YZ Pressure Map');
% 
% % pcolor with shading interp is excellent for 2D smooth heatmaps
% pcolor(Y_grid_fine, Z_grid_fine, P_fine);
% shading interp; % Removes grid lines and blends colors smoothly
% 
% % Formatting the plot
% colormap(turbo); % 'turbo' or 'jet' are great for acoustic pressure maps
% cb = colorbar;
% ylabel(cb, 'Pressure (Arbitrary Units)', 'FontSize', 11, 'FontWeight', 'bold');
% 
% title('PMUT Beam Shape: Smoothed YZ Pressure Map', 'FontSize', 14);
% xlabel('Lateral Position Y (mm)', 'FontSize', 12);
% ylabel('Depth Z (mm)', 'FontSize', 12);
% 
% % Reverse the Z-axis (plotted on Y) so depth increases as you go down the chart
% set(gca, 'YDir', 'reverse', 'FontSize', 11); 
% axis image;



% ==========================================
% --- CONFIGURATION ---
% ==========================================
% filename = 'C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\YZ_focused_detailed\Scan_2D_UMSmap.txt';
filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\XZ_detailed\Scan_2D_UMSmap.txt";
scanPlane = 'XZ';                % Set to 'XZ' or 'YZ' to update plot labels
showFocalZone = false;            % Set to true to display the -6dB contour, false to hide
% ==========================================

% 1. Load the Data
% The file has 4 header lines before the numeric data
rawData = readmatrix(filename, 'NumHeaderLines', 4);

% 2. Extract Axes and Pressure Matrix dynamically
% This automatically adjusts to however many columns your file has (16 for XZ, 12 for YZ)
Lat_raw = rawData(1, 2:end);  % Lateral coordinates
Z_raw = rawData(2:end, 1);    % Depth coordinates
P_raw = rawData(2:end, 2:end); % Pressure matrix

% 3. Coordinate Adjustments
% Center the lateral axis around 0 mm
Lat = Lat_raw - mean(Lat_raw);

% Shift the Z-axis so the scan starts exactly at 3 mm depth
Z = Z_raw - Z_raw(1) + 3;

% 4. Interpolation for Smoothing
[Lat_grid, Z_grid] = meshgrid(Lat, Z);

% Create a fine grid for smooth visualization
Lat_fine = linspace(min(Lat), max(Lat), 200);
Z_fine = linspace(min(Z), max(Z), 200);
[Lat_grid_fine, Z_grid_fine] = meshgrid(Lat_fine, Z_fine);

% Smooth the pressure data
P_fine = interp2(Lat_grid, Z_grid, P_raw, Lat_grid_fine, Z_grid_fine, 'spline');

% Calculate Decibel (dB) map for the focal zone
% Acoustic pressure uses 20*log10(P / P_max)
P_dB = 20 * log10(P_fine ./ max(P_fine(:)));
P_lin = P_fine ./ max(P_fine(:));

% 5. Plotting the Heatmap
figure('Color', 'w', 'Name', ['PMUT ', scanPlane, ' Pressure Map']);

pcolor(Lat_grid_fine, Z_grid_fine, P_fine);
shading interp;
colormap(turbo); 
cb = colorbar;
ylabel(cb, 'Pressure (Arbitrary Units)', 'FontSize', 11, 'FontWeight', 'bold');

% 6. Toggle -6dB Focal Zone
if showFocalZone
    hold on;
    % Plot the -6dB contour line in a dashed white line so it stands out
    [~, h] = contour(Lat_grid_fine, Z_grid_fine, P_dB, [-6 -6], 'w--', 'LineWidth', 2);
    
    % Add a custom legend just for the focal zone
    legend(h, '-6dB Focal Zone', 'Location', 'northeast', 'TextColor', 'w', ...
        'Color', 'none', 'EdgeColor', 'none', 'FontSize', 10);
    hold off;
end

% 7. Dynamic Formatting
title(['PMUT Beam Shape: Smoothed ', scanPlane, ' Pressure Map'], 'FontSize', 14);
% Automatically uses 'X' or 'Y' for the label based on your scanPlane setting
xlabel(['Lateral Position ', scanPlane(1), ' (mm)'], 'FontSize', 12); 
ylabel('Depth Z (mm)', 'FontSize', 12);

% Reverse Z-axis so depth goes downwards
set(gca, 'YDir', 'reverse', 'FontSize', 11); 
axis tight; % Keeps the axes completely snug against your data borders
pbaspect([1 2 1]); % Forces the plot box to a specific [Width Height 1] ratio


% ---------------------------- db map ---------------------------------

figure('Color', 'w', 'Name', ['PMUT ', scanPlane, ' Pressure Map']);

pcolor(Lat_grid_fine, Z_grid_fine, P_dB);
shading interp;
colormap(turbo); 
cb = colorbar;
ylabel(cb, 'Pressure (Arbitrary Units)', 'FontSize', 11, 'FontWeight', 'bold');

title(['PMUT Beam Shape: Smoothed ', scanPlane, ' Pressure Map'], 'FontSize', 14);
% Automatically uses 'X' or 'Y' for the label based on your scanPlane setting
xlabel(['Lateral Position ', scanPlane(1), ' (mm)'], 'FontSize', 12); 
ylabel('Depth Z (mm)', 'FontSize', 12);

% Reverse Z-axis so depth goes downwards
set(gca, 'YDir', 'reverse', 'FontSize', 11); 
axis tight; % Keeps the axes completely snug against your data borders
pbaspect([1 2 1]); % Forces the plot box to a specific [Width Height 1] ratio