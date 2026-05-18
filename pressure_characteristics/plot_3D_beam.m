% % ==========================================
% % --- CONFIGURATION ---
% % ==========================================
% filename_XZ = 'C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\YZ_focused_detailed\Scan_2D_UMSmap.txt'; % Replace with your actual XZ file name
% filename_YZ = "C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\XZ_detailed\Scan_2D_UMSmap.txt"; % Replace with your actual YZ file name
% noiseThreshold = 0.65;       % Any pressure below 5% of max is considered "empty space"
% beamPointSize = 25;          % Uniform size for the main beam points
% % ==========================================
% 
% % 1. Load the Data
% data_XZ = readmatrix(filename_XZ, 'NumHeaderLines', 4);
% data_YZ = readmatrix(filename_YZ, 'NumHeaderLines', 4);
% 
% X = data_XZ(1, 2:end) - mean(data_XZ(1, 2:end));
% Z_XZ = data_XZ(2:end, 1) - data_XZ(2, 1) + 3;
% P_XZ = data_XZ(2:end, 2:end);
% 
% Y = data_YZ(1, 2:end) - mean(data_YZ(1, 2:end));
% Z_YZ = data_YZ(2:end, 1) - data_YZ(2, 1) + 3;
% P_YZ = data_YZ(2:end, 2:end);
% 
% % 2. Create Common Fine Grids
% % We use a 60x60x100 resolution. Going too high (e.g., 200x200x200) 
% % will create 8 million points and freeze your computer!
% X_fine = linspace(min(X), max(X), 60);
% Y_fine = linspace(min(Y), max(Y), 60);
% Z_fine = linspace(min(min(Z_XZ), min(Z_YZ)), max(max(Z_XZ), max(Z_YZ)), 100);
% 
% % 3. Interpolate the 2D planes to match the common Z-axis
% [X_grid, Z_XZ_grid] = meshgrid(X, Z_XZ);
% [X_fine_grid, Z_fine_grid_X] = meshgrid(X_fine, Z_fine);
% P_XZ_interp = interp2(X_grid, Z_XZ_grid, P_XZ, X_fine_grid, Z_fine_grid_X, 'spline');
% 
% [Y_grid, Z_YZ_grid] = meshgrid(Y, Z_YZ);
% [Y_fine_grid, Z_fine_grid_Y] = meshgrid(Y_fine, Z_fine);
% P_YZ_interp = interp2(Y_grid, Z_YZ_grid, P_YZ, Y_fine_grid, Z_fine_grid_Y, 'spline');
% 
% % Ensure no negative pressures from spline overshoot
% P_XZ_interp(P_XZ_interp < 0) = 0;
% P_YZ_interp(P_YZ_interp < 0) = 0;
% 
% % 4. Construct the Full 3D Volume
% [X_3D, Y_3D, Z_3D] = meshgrid(X_fine, Y_fine, Z_fine);
% P_3D = zeros(size(X_3D));
% 
% % Multiply the X and Y profiles at every depth (Z) to form the 3D volume
% for k = 1:length(Z_fine)
%     prof_X = P_XZ_interp(k, :);
%     prof_Y = P_YZ_interp(k, :);
%     % Outer product projects the 1D profiles into a 2D horizontal slice
%     P_3D(:, :, k) = (prof_Y') * prof_X; 
% end
% 
% % Normalize the whole 3D volume from 0 to 1
% P_3D = P_3D ./ max(P_3D(:));
% 
% % Flatten arrays for the scatter plot
% X_vec = X_3D(:); Y_vec = Y_3D(:); Z_vec = Z_3D(:); P_vec = P_3D(:);
% 
% % 5. Separate "Beam" from "Empty Space"
% isBeam = P_vec >= noiseThreshold;
% isBg = P_vec < noiseThreshold;
% 
% % 6. Plotting
% figure('Color', 'w', 'Name', 'Full 3D Beam Volume');
% 
% % First, plot the empty background space as tiny, faint grey dots
% scatter3(X_vec(isBg), Y_vec(isBg), Z_vec(isBg), 1, [0.8 0.8 0.8], ...
%     'filled', 'MarkerFaceAlpha', 0.05);
% hold on;
% 
% % Second, plot the actual beam points: uniform size, colored by pressure
% scatter3(X_vec(isBeam), Y_vec(isBeam), Z_vec(isBeam), beamPointSize, P_vec(isBeam), ...
%     'filled', 'MarkerFaceAlpha', 0.4); % Alpha 0.4 lets you see through the outer layers
% hold off;
% 
% % 7. Formatting
% colormap(turbo);
% cb = colorbar;
% ylabel(cb, 'Normalized Pressure', 'FontSize', 11, 'FontWeight', 'bold');
% 
% title('PMUT 3D Volumetric Beam Reconstruction', 'FontSize', 14);
% xlabel('X Position (mm)', 'FontSize', 12);
% ylabel('Y Position (mm)', 'FontSize', 12);
% zlabel('Depth Z (mm)', 'FontSize', 12);
% 
% set(gca, 'ZDir', 'reverse', 'FontSize', 11);
% view([-35, 30]); 
% grid on;
% pbaspect([1 1 2]);


%-------------------------------- V2 --------------------------------------

% % ==========================================
% % --- CONFIGURATION ---
% % ==========================================
% filename_XZ = 'C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\YZ_focused_detailed\Scan_2D_UMSmap.txt'; % Replace with your actual XZ file name
% filename_YZ = "C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\XZ_detailed\Scan_2D_UMSmap.txt"; % Replace with your actual YZ file name
% noiseThreshold = 0.60;       % Removes points below 5% pressure to keep the plot clean
% % ==========================================
% 
% % 1. Load the Data
% data_XZ = readmatrix(filename_XZ, 'NumHeaderLines', 4);
% data_YZ = readmatrix(filename_YZ, 'NumHeaderLines', 4);
% 
% X = data_XZ(1, 2:end) - mean(data_XZ(1, 2:end));
% Z_XZ = data_XZ(2:end, 1) - data_XZ(2, 1) + 3;
% P_XZ = data_XZ(2:end, 2:end);
% 
% Y = data_YZ(1, 2:end) - mean(data_YZ(1, 2:end));
% Z_YZ = data_YZ(2:end, 1) - data_YZ(2, 1) + 3;
% P_YZ = data_YZ(2:end, 2:end);
% 
% % 2. Create Common Fine Grids
% X_fine = linspace(min(X), max(X), 60);
% Y_fine = linspace(min(Y), max(Y), 60);
% Z_fine = linspace(min(min(Z_XZ), min(Z_YZ)), max(max(Z_XZ), max(Z_YZ)), 100);
% 
% % 3. Interpolate the 2D planes to match the common Z-axis
% [X_grid, Z_XZ_grid] = meshgrid(X, Z_XZ);
% [X_fine_grid, Z_fine_grid_X] = meshgrid(X_fine, Z_fine);
% P_XZ_interp = interp2(X_grid, Z_XZ_grid, P_XZ, X_fine_grid, Z_fine_grid_X, 'spline');
% 
% [Y_grid, Z_YZ_grid] = meshgrid(Y, Z_YZ);
% [Y_fine_grid, Z_fine_grid_Y] = meshgrid(Y_fine, Z_fine);
% P_YZ_interp = interp2(Y_grid, Z_YZ_grid, P_YZ, Y_fine_grid, Z_fine_grid_Y, 'spline');
% 
% % Ensure no negative pressures from spline overshoot
% P_XZ_interp(P_XZ_interp < 0) = 0;
% P_YZ_interp(P_YZ_interp < 0) = 0;
% 
% % 4. Construct the Full 3D Volume
% [X_3D, Y_3D, Z_3D] = meshgrid(X_fine, Y_fine, Z_fine);
% P_3D = zeros(size(X_3D));
% 
% for k = 1:length(Z_fine)
%     prof_X = P_XZ_interp(k, :);
%     prof_Y = P_YZ_interp(k, :);
%     P_3D(:, :, k) = (prof_Y') * prof_X; 
% end
% 
% % Normalize the whole 3D volume from 0 to 1
% P_3D = P_3D ./ max(P_3D(:));
% 
% % Flatten arrays for the scatter plot
% X_vec = X_3D(:); Y_vec = Y_3D(:); Z_vec = Z_3D(:); P_vec = P_3D(:);
% 
% % 5. Filter out the absolute lowest noise for a clean paper figure
% isValid = P_vec >= noiseThreshold;
% X_plot = X_vec(isValid);
% Y_plot = Y_vec(isValid);
% Z_plot = Z_vec(isValid);
% P_plot = P_vec(isValid);
% 
% % 6. Calculate Dynamic Point Sizes
% minPointSize = 1;
% maxPointSize = 120; % Adjust this if the core looks too blocky or too thin
% 
% % Cubing the pressure (P_plot .^ 3) forces low values to stay very small 
% % and only the highest values to balloon up in size.
% pointSizes = minPointSize + (maxPointSize - minPointSize) * (P_plot .^ 3);
% 
% % 7. Plotting
% figure('Color', 'w', 'Name', 'Full 3D Beam Volume');
% 
% scatter3(X_plot, Y_plot, Z_plot, pointSizes, P_plot, 'filled', 'MarkerFaceAlpha', 0.4);
% 
% % 8. Formatting
% colormap(turbo);
% cb = colorbar;
% ylabel(cb, 'Normalized Pressure', 'FontSize', 11, 'FontWeight', 'bold');
% 
% title('PMUT 3D Volumetric Beam Reconstruction', 'FontSize', 14);
% xlabel('X Position (mm)', 'FontSize', 12);
% ylabel('Y Position (mm)', 'FontSize', 12);
% zlabel('Depth Z (mm)', 'FontSize', 12);
% 
% set(gca, 'ZDir', 'reverse', 'FontSize', 11);
% view([-35, 30]); 
% grid on;
% pbaspect([1 1 2]);

% -------------------------------- V3 -------------------------------------

% ==========================================
% --- CONFIGURATION ---
% ==========================================
filename_XZ = 'C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\YZ_focused_detailed\Scan_2D_UMSmap.txt'; % Replace with your actual XZ file name
filename_YZ = "C:\Users\cioca100\Desktop\WaterTank_Imaging\Array27\XZ_detailed\Scan_2D_UMSmap.txt"; % Replace with your actual YZ file name
noiseThreshold = 0.65;        % Try bumping this to 0.4 or 0.5 to isolate the core!
% ==========================================

% 1. Load the Data
data_XZ = readmatrix(filename_XZ, 'NumHeaderLines', 4);
data_YZ = readmatrix(filename_YZ, 'NumHeaderLines', 4);

X = data_XZ(1, 2:end) - mean(data_XZ(1, 2:end));
Z_XZ = data_XZ(2:end, 1) - data_XZ(2, 1) + 3;
P_XZ = data_XZ(2:end, 2:end);

Y = data_YZ(1, 2:end) - mean(data_YZ(1, 2:end));
Z_YZ = data_YZ(2:end, 1) - data_YZ(2, 1) + 3;
P_YZ = data_YZ(2:end, 2:end);

% 2. Create Common Fine Grids
X_fine = linspace(min(X), max(X), 60);
Y_fine = linspace(min(Y), max(Y), 60);
Z_fine = linspace(min(min(Z_XZ), min(Z_YZ)), max(max(Z_XZ), max(Z_YZ)), 100);

% 3. Interpolate the 2D planes to match the common Z-axis
[X_grid, Z_XZ_grid] = meshgrid(X, Z_XZ);
[X_fine_grid, Z_fine_grid_X] = meshgrid(X_fine, Z_fine);
P_XZ_interp = interp2(X_grid, Z_XZ_grid, P_XZ, X_fine_grid, Z_fine_grid_X, 'spline');

[Y_grid, Z_YZ_grid] = meshgrid(Y, Z_YZ);
[Y_fine_grid, Z_fine_grid_Y] = meshgrid(Y_fine, Z_fine);
P_YZ_interp = interp2(Y_grid, Z_YZ_grid, P_YZ, Y_fine_grid, Z_fine_grid_Y, 'spline');

P_XZ_interp(P_XZ_interp < 0) = 0;
P_YZ_interp(P_YZ_interp < 0) = 0;

% 4. Construct the Full 3D Volume
[X_3D, Y_3D, Z_3D] = meshgrid(X_fine, Y_fine, Z_fine);
P_3D = zeros(size(X_3D));

for k = 1:length(Z_fine)
    prof_X = P_XZ_interp(k, :);
    prof_Y = P_YZ_interp(k, :);
    P_3D(:, :, k) = (prof_Y') * prof_X; 
end

% Normalize the whole 3D volume from 0 to 1
P_3D = P_3D ./ max(P_3D(:));

% Flatten arrays for the scatter plot
X_vec = X_3D(:); Y_vec = Y_3D(:); Z_vec = Z_3D(:); P_vec = P_3D(:);

% 5. Filter out points below the threshold
isValid = P_vec >= noiseThreshold;
X_plot = X_vec(isValid);
Y_plot = Y_vec(isValid);
Z_plot = Z_vec(isValid);
P_plot = P_vec(isValid);

% 6. Calculate Dynamic Point Sizes
minPointSize = 1;
maxPointSize = 120; 

pointSizes = minPointSize + (maxPointSize - minPointSize) * (P_plot .^ 3);

% 7. Plotting
figure('Color', 'w', 'Name', 'Full 3D Beam Volume Isolated');

scatter3(X_plot, Y_plot, Z_plot, pointSizes, P_plot, 'filled', 'MarkerFaceAlpha', 0.4);

% 8. Formatting & Color Locking
colormap(turbo);

% --- THE FIX: Hard-lock the color scale to absolute normalized pressure ---
clim([0 1]); % Note: If you are using MATLAB R2021b or older, change 'clim' to 'caxis'
% --------------------------------------------------------------------------

cb = colorbar;
ylabel(cb, 'Normalized Pressure', 'FontSize', 11, 'FontWeight', 'bold');

title(sprintf('PMUT 3D Volumetric Beam (Threshold: >%.2f)', noiseThreshold), 'FontSize', 14);
xlabel('X Position (mm)', 'FontSize', 12);
ylabel('Y Position (mm)', 'FontSize', 12);
zlabel('Depth Z (mm)', 'FontSize', 12);

set(gca, 'ZDir', 'reverse', 'FontSize', 11);
view([-35, 30]); 
grid on;
pbaspect([1 1 3]);



% --- Assuming X_3D, Y_3D, Z_3D, and P_3D are already computed ---

% Convert the normalized pressure to Decibels (dB) for standard acoustic mapping
P_3D_dB = 20 * log10(P_3D);

% Create a new, clean figure for the paper
figure('Color', 'w', 'Name', 'PMUT 3D Isosurface Beam', 'Position', [100, 100, 800, 600]);

hold on;

% 1. Plot the Outer Shell (-6dB boundary)
% This represents the standard width of the beam
[faces_out, verts_out] = isosurface(X_3D, Y_3D, Z_3D, P_3D_dB, -6);
p1 = patch('Vertices', verts_out, 'Faces', faces_out);
p1.FaceColor = [0.2 0.6 1.0]; % A nice cool blue
p1.EdgeColor = 'none';        % Remove jagged mesh lines
p1.FaceAlpha = 0.3;           % Make it 30% transparent so we can see inside!

% 2. Plot the Inner Focal Core (-3dB boundary)
% This represents the highest intensity region
[faces_in, verts_in] = isosurface(X_3D, Y_3D, Z_3D, P_3D_dB, -3);
p2 = patch('Vertices', verts_in, 'Faces', faces_in);
p2.FaceColor = [1.0 0.3 0.2]; % A bold warm red/orange
p2.EdgeColor = 'none';
p2.FaceAlpha = 0.9;           % Make it almost completely solid

hold off;

% 3. Add Lighting and Shadows (Crucial for 2D Snapshots!)
camlight('headlight');   % Shines a light from the camera's perspective
camlight('left');        % Adds a secondary fill light from the side
lighting gouraud;        % Calculates smooth shadows across the curved surfaces
material dull;           % Prevents the surface from looking too shiny/plastic

% 4. Formatting the Plot
title('PMUT 3D Volumetric Beam Profile', 'FontSize', 14);
xlabel('X Position (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y Position (mm)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Depth Z (mm)', 'FontSize', 12, 'FontWeight', 'bold');

% Standardize the view
set(gca, 'ZDir', 'reverse', 'FontSize', 11);
view([-35, 25]); % A slightly lower angle usually looks better for isosurfaces
grid on;
box on; % Puts a bounding box around the 3D space to anchor it visually
pbaspect([1 1 2]); 

% Add a legend to explain the shells
legend([p1, p2], {'-6dB Boundary', '-3dB Focal Core'}, 'Location', 'northeast');