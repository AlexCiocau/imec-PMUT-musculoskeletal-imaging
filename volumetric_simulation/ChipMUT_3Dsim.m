%% Volumetric Huygens-Fresnel simulation of a channel-addressable 2D PMUT array
%
% Models every cell as a coherent point source at its real (x, y, 0) position.
% Applies a single transmit delay per CHANNEL (all cells in a column share it),
% which is what causes the array to be astigmatic: azimuth (x) is electronically
% focused but elevation (y) only has its natural focus z_NF,y = D_y^2 / (4 lambda).
%
% Two visualisations:
%   (a) 3D scatter of every voxel whose intensity is above THRESHOLD_DB.
%   (b) Azimuth-axial and elevation-axial slices through the on-axis line.
%
% No toolboxes required.

clear; clc; close all;

%% ---------- 1. Acoustic parameters ------------------------------------
c       = 1480;             % speed of sound in water (m/s)
fc      = 10.5e6;            % center frequency (Hz)
lambda  = c / fc;
k       = 2*pi / lambda;

%% ---------- 2. Array geometry -----------------------------------------
% Channels are arranged along x (azimuth). Each channel is a column of cells
% along y (elevation). If your real cell layout is different (e.g. 8 x 10
% instead of 1 x 80), change N_cell_per_ch_x and N_cell_per_ch_y.
N_ch              = 64;
ch_pitch          = 75e-6;        % channel pitch (azimuth)
N_cell_per_ch_x   = 1;            % cells per channel along x: 64
N_cell_per_ch_y   = 80;           % cells per channel along y: 80
cell_pitch_x      = 75e-6;
cell_pitch_y      = 75e-6;        % <-- change this to match your real chip

% Build the (x, y) position of every cell, plus its parent channel index
[ix, iy, ich]     = ndgrid(1:N_cell_per_ch_x, ...
                            1:N_cell_per_ch_y, ...
                            1:N_ch);
ich               = ich(:);
ch_centre_x       = (ich - (N_ch+1)/2) * ch_pitch;
dx_in_ch          = (ix(:) - (N_cell_per_ch_x+1)/2) * cell_pitch_x;
dy_in_ch          = (iy(:) - (N_cell_per_ch_y+1)/2) * cell_pitch_y;

x_cell            = ch_centre_x + dx_in_ch;
y_cell            = dy_in_ch;
N_cells           = numel(x_cell);

% Apertures and natural focal depths
D_x  = (N_ch-1)*ch_pitch + N_cell_per_ch_x*cell_pitch_x;
D_y  = (N_cell_per_ch_y-1) * cell_pitch_y;
zNFx = D_x^2 / (4*lambda);
zNFy = D_y^2 / (4*lambda);

fprintf('Total cells       : %d (%d ch x %d cells/ch)\n', ...
        N_cells, N_ch, N_cell_per_ch_x*N_cell_per_ch_y);
fprintf('Aperture (Dx, Dy) : %.2f mm  x  %.2f mm\n', D_x*1e3, D_y*1e3);
fprintf('Natural focus     : zNF,x = %.1f mm,  zNF,y = %.1f mm\n', ...
        zNFx*1e3, zNFy*1e3);

%% ---------- 3. Choose transmit delay law ------------------------------
% NOTE: focal gain saturates for z_focus > zNF,x. If you set z_focus past
% the natural focus, the on-axis intensity peak will *not* be at z_focus -
% it will sit near zNF,x. The simulation does this correctly because it is
% pure physics; nothing in it knows that you "asked" for a deeper focus.
mode    = 'plane';      % 'plane' or 'focus'
x_focus = 0;
z_focus = 25e-3;        % try 20e-3, 50e-3, 70e-3 to see focal-gain saturation

ch_x_centres = ((1:N_ch) - (N_ch+1)/2) * ch_pitch;
switch mode
    case 'plane'
        ch_tau = zeros(1, N_ch);
    case 'focus'
        d_ch    = sqrt((x_focus - ch_x_centres).^2 + z_focus^2);
        ch_tau  = (max(d_ch) - d_ch) / c;
    otherwise
        error('mode must be ''plane'' or ''focus''.');
end
cell_tau = ch_tau(ich);                     % cells inherit channel delay

%% ---------- 4. Imaging volume ----------------------------------------
% Resolution-vs-runtime knobs. With these defaults and 5120 cells the
% computation takes ~30-60 s in regular MATLAB. Halve N_z first if too slow.
Nx = 50;  Ny = 50;  Nz = 120;
x_grid = single(linspace(-6e-3,  6e-3, Nx));
y_grid = single(linspace(-6e-3,  6e-3, Ny));
% z_grid = single(linspace( 0.5e-3, 80e-3, Nz));
z_grid = single(linspace( 0.5e-3, 120e-3, Nz));
[X3, Y3, Z3] = ndgrid(x_grid, y_grid, z_grid);

%% ---------- 5. Coherent sum over all cells ----------------------------
% Use single precision for speed and memory. Each iteration adds one cell
% to the running complex pressure P3.
x_cell_s   = single(x_cell);
y_cell_s   = single(y_cell);
cell_tau_s = single(cell_tau);
k_s        = single(k);
omega_s    = single(2*pi*fc);

P3 = complex(zeros(Nx, Ny, Nz, 'single'));
fprintf('Computing 3D field over %dx%dx%d voxels and %d cells ...\n', ...
        Nx, Ny, Nz, N_cells);
tic;
for i = 1:N_cells
    if mod(i, 512) == 0
        fprintf('  cell %5d / %d   (%.1f s elapsed)\n', i, N_cells, toc);
    end
    r     = sqrt((X3 - x_cell_s(i)).^2 + (Y3 - y_cell_s(i)).^2 + Z3.^2);
    phase = -(k_s*r + omega_s*cell_tau_s(i));
    P3    = P3 + exp(1i*phase) ./ r;
end
fprintf('Done (%.1f s).\n', toc);

I3      = abs(P3).^2;
I3_dB   = 10*log10(I3 / max(I3(:)) + eps);

%% ---------- 6. 3D scatter above intensity threshold -------------------
% Only voxels whose intensity is above THRESHOLD_DB (relative to peak) are
% drawn. Lower the threshold to see more of the side-lobe structure;
% raise it to isolate the bright focal core.
THRESHOLD_DB = -10;

mask = I3_dB > THRESHOLD_DB;
xv = double(X3(mask))*1e3;
yv = double(Y3(mask))*1e3;
zv = double(Z3(mask))*1e3;
cv = double(I3_dB(mask));
fprintf('Voxels above %.1f dB: %d / %d (%.2f %%)\n', ...
        THRESHOLD_DB, nnz(mask), numel(mask), 100*nnz(mask)/numel(mask));

figure('Color','w','Position',[80 80 820 620]);
scatter3(xv, yv, zv, 14, cv, 'filled', 'MarkerEdgeAlpha', 0);
clim([THRESHOLD_DB 0]);
colormap(turbo);
cb = colorbar; cb.Label.String = 'Intensity (dB)';
xlabel('Azimuth x (mm)');
ylabel('Elevation y (mm)');
zlabel('Depth z (mm)');
set(gca, 'ZDir','reverse');
% axis equal vis3d; 
grid on;
view(135, 25);
pbaspect([1 1 2]); % Forces the plot box to a specific [Width Height 1] ratio
title(sprintf('Voxels above %.0f dB  (mode = %s)', THRESHOLD_DB, mode));

%% ---------- 7. Azimuth-axial and elevation-axial slices ---------------
[~, iy0]   = min(abs(y_grid));
[~, ix0]   = min(abs(x_grid));
slice_xz   = squeeze(I3_dB(:, iy0, :));
slice_yz   = squeeze(I3_dB(ix0, :, :));

figure('Color','w','Position',[920 80 720 720]);
subplot(2,1,1);
imagesc(double(x_grid)*1e3, double(z_grid)*1e3, double(slice_xz')); clim([-20 0]);
% colormap(hot); 
colorbar;
xlabel('Azimuth x (mm)'); ylabel('Depth z (mm)');
title('Azimuth-axial plane (y = 0) - electronically focused direction');
axis image; hold on;
yline(zNFx*1e3, 'g--', 'LineWidth', 1, 'Label', 'z_{NF,x}', 'Color', [0 .8 0]);
if strcmp(mode,'focus')
    plot(x_focus*1e3, z_focus*1e3, 'g+', 'MarkerSize', 14, 'LineWidth', 2);
end
hold off;

subplot(2,1,2);
imagesc(double(y_grid)*1e3, double(z_grid)*1e3, double(slice_yz')); clim([-20 0]);
colormap(hot); colorbar;
xlabel('Elevation y (mm)'); ylabel('Depth z (mm)');
title(sprintf('Elevation-axial plane (x = 0) - natural focus only at z_{NF,y} = %.1f mm', zNFy*1e3));
axis image; hold on;
yline(zNFy*1e3, 'g--', 'LineWidth', 1, 'Label', 'z_{NF,y}', 'Color', [0 .8 0]);
hold off;

%% ---------- 8. On-axis intensity in both lateral directions ----------
% Useful for measuring the astigmatic separation: the depth at which the
% azimuth slice peaks vs the depth at which the elevation slice peaks.
on_axis_x = squeeze(I3_dB(ix0, iy0, :));
[~, idx_max] = max(on_axis_x);
fprintf('On-axis intensity peaks at z = %.1f mm (you asked for z_focus = %.1f mm)\n', ...
        z_grid(idx_max)*1e3, z_focus*1e3);

figure('Color','w','Position',[80 720 720 320]);
plot(double(z_grid)*1e3, double(on_axis_x), 'LineWidth', 1.4); grid on;
xlabel('Depth z (mm)');
ylabel('On-axis intensity (dB)');
title('On-axis intensity vs depth');
hold on;
xline(zNFx*1e3, 'g--', 'LineWidth', 1.2);
xline(zNFy*1e3, 'b--', 'LineWidth', 1.2);
if strcmp(mode,'focus')
    xline(z_focus*1e3, 'r--', 'LineWidth', 1.2);
end
legend('on-axis intensity', 'z_{NF,x}', 'z_{NF,y}', 'requested z_{focus}', ...
       'Location','best');
hold off;