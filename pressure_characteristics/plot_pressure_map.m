% % Define the filename (make sure the file is in your current MATLAB directory)
% filename = 'Scan_2D_UMSmap.txt';
% 
% % Read the numeric data, skipping the first 4 header lines
% % The readmatrix function will handle the tab-delimited format automatically
% data = readmatrix(filename, 'NumHeaderLines', 4);

% 
% % Extract X and Y axes
% % The first row (cols 2 to end) contains the X coordinates
% X = data(1, 2:end);
% 
% % The first column (rows 2 to end) contains the Y coordinates
% Y = data(2:end, 1);
% 
% % Extract the 11x11 measurement matrix
% Z = data(2:end, 2:end);
% 
% % Normalize the Z data (Min-Max normalization to a 0 to 1 scale)
% Z_min = min(Z(:));
% Z_max = max(Z(:));
% Z_norm = (Z - Z_min) / (Z_max - Z_min);
% 
% % Create the figure and plot the normalized heatmap
% figure;
% % imagesc is great for plotting gridded matrix data
% imagesc(X, Y, Z_norm);
% 
% % Set Y-axis to increase upwards (standard for physical spatial mapping)
% set(gca, 'YDir', 'normal'); 
% 
% % Add a colorbar to show the normalized scale
% colorbar;
% % You can change 'parula' to 'jet', 'hot', or 'turbo' depending on your preference
% colormap('turbo'); 
% 
% % Add labels and title
% xlabel('X Position');
% ylabel('Y Position');
% title('Normalized PMUT Array Surface Pressure Map');








% Define the filename 
filename = 'Scan_2D_UMSmap.txt';

% Apparently Array27 ?uneven? Comments: "Array27 - SurfaceMap @10.5MHz,
% depth 3.245mm, aaligned visually; XY - 11x11, res 0.7x0.6"
filename = "C:\Users\cioca100\Desktop\Final_Measurements\Array27\SurfacePressureMap\Scan_2D_UMSmap.txt";

% Array 20 Surface Pressure
% filename  = "C:\Users\cioca100\Desktop\WaterTank_Imaging\Array20\SurfaceMap\Scan_2D_UMSmap.txt";

% Array 26 Surface Pressure
% filename  = "C:\Users\cioca100\Desktop\WaterTank_Imaging\Array26\better_SurfaceMap\Scan_2D_UMSmap.txt";

% Array 27 Surface Pressure
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\Surface_pressure\Scan_2D_UMSmap.txt";

% Array 47 Surface Pressure
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\Array47\SurfaceMap\Scan_2D_UMSmap.txt";
% Read the numeric data, skipping the first 4 header lines
data = readmatrix(filename, 'NumHeaderLines', 4);

% Extract X and Y axes
X = data(1, 2:end);
Y = data(2:end, 1);

% Extract the 11x11 measurement matrix
Z = data(2:end, 2:end);

% Normalize the Z data (Min-Max normalization to a 0 to 1 scale)
Z_min = min(Z(:));
Z_max = max(Z(:));
Z_norm = (Z - Z_min) / (Z_max - Z_min);

% --- INTERPOLATION SETUP ---

% Create a grid for the original 11x11 data points
[X_grid, Y_grid] = meshgrid(X, Y);

% Create a finer grid for a smoother high-resolution plot (e.g., 100x100 points)
X_fine = linspace(min(X), max(X), 80);
Y_fine = linspace(min(Y), max(Y), 64);
[X_fine_grid, Y_fine_grid] = meshgrid(X_fine, Y_fine);

% Interpolate the normalized Z data over the finer grid using cubic spline interpolation
Z_interp = interp2(X_grid, Y_grid, Z_norm, X_fine_grid, Y_fine_grid, 'spline');

% --- PLOTTING ---

figure;
% Plot the smoothly interpolated data
imagesc(X_fine, Y_fine, Z_interp);

% Set Y-axis to increase upwards
set(gca, 'YDir', 'normal'); 

% Add a colorbar and choose a good colormap for heatmaps
colorbar;
colormap('turbo'); % 'turbo' provides great contrast for pressure maps, but 'parula' or 'jet' work too!

% Add labels and title
xlabel('X Position');
ylabel('Y Position');
title('Smooth Normalized PMUT Array Surface Pressure Map');
axis image;