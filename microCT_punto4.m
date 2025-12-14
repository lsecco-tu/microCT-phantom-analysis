%% AP3232 - MEDICAL IMAGING - ASSIGNMENT
%  BEAM HARDENING IN MICRO CT
%  Micro CT scanner experiment 
%
%  This code allows for the visualization and analysis of the reconstructed volume of the
%  phantom under consideration
%

% Save results
outDir = uigetdir(pwd, 'Select folder to save Command Window output');

if outDir == 0
    disp('User cancelled folder selection');
    return
end

% Define output file
logFile = fullfile(outDir, 'command_window_log.txt');

% Start recording Command Window
diary(logFile);
diary on

% Data
file_name = 'CT_24_01_01.nii';
data = load_nii(file_name);
img_uint16 = data.img;                           % 3D array (class uint16)
[dim_x, dim_y, dim_z] = size(img_uint16);        % sizes of 3D array

% Conversion of voxel values to attenuation coefficients
% Range of attenuation coefficients
mu_min = -0.0002;   % [mm^-1]
mu_max = 0.02;      % [mm^-1]

img_double = double(img_uint16);

% Linear mapping 0..65535 -> mu_min...mu_max
mu_volume = mu_min + (img_double / 65535) * (mu_max - mu_min);

% Correction of reconstruction sotware bug
mu_volume_pi = mu_volume * (2*pi);

% figure setup
set(groot, 'defaultAxesFontSize', 18);
set(groot, 'defaultTextFontSize', 24);


%% 1. VISUALIZATION OF PHANTOM SECTIONS IN THE TRANSVERSAL PLANE (x,y)

% Slice indexes increase from the bottom of the phantom (z = 1, 55 mm
% plexi), to the top of the phantom (z = dim_z, 15 mm plexi)

if strcmp(file_name, 'CT_21_1_01.nii')
    slices = [350, 250, 150, 80, 20];   % 15 mm, 35 mm, 55 mm, 55 mm (2 rods), 55 mm (3 rods)
elseif strcmp(file_name, 'CT_24_01_01.nii')
    slices = [420, 350, 250, 180, 90];
end

diameters = [15, 35, 55];
n_slices = 30;                             % number of slices used to reduce noise
all_avg_images = cell(1, length(slices));

for i = 1:length(slices)

    slice = slices(i);        % current transversal section

    % Volume of slices to average to reduce noise
    slices_to_average = mu_volume_pi(:, :, slice:(slice + n_slices));
    avg_image = mean(slices_to_average,3);

    % Averaged transversal section images stored in a cell array
    all_avg_images{i} = avg_image;

    % Display of current averaged image
    figure(); imagesc(avg_image); axis equal; colormap gray; colorbar; axis image;
    title(['Transversal Plane (Z = ' num2str(slice) ')']);
    xlabel('Y - vertical direction'); 
    ylabel('X - horizontal direction');  

end

%% 2. LINE PROFILES THROUGH CENTER (plexi-only sections)

centers = zeros(3,2); % array the stores the center coordinates for each section 
                      % (in order: 15 mm, 35 mm, 55 mm)

for i = 1:3   % indexes of plexi-only sections (slice_indexes array)

    current_avg = all_avg_images{i};   % current averaged image

    % Calculation of plexi center coordinates
    % threshold to separate section from background
    t = (min(current_avg(:)) + max(current_avg(:))) / 2; 
    mask = current_avg > t;
    [cy, cx] = circleCenterMask(mask);      %centre of plexi section
                                            %y is the vertical direction
                                            %(cy is the x coordinate in the
                                            %image)
    centers(i,:) = [round(cy), round(cx)]; 

    line_profile = current_avg(:, round(cy));   %centre row → along y direction

    % Visualize section, center and line through center
    figure; imagesc(current_avg); axis image; colorbar; colormap gray; hold on;
    plot(cy, cx, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    plot([cy cy], [1 size(current_avg,1)], 'r--', 'LineWidth', 1);
    title('Line through calculated center');
    xlabel('Y - vertical direction')
    ylabel('X - horizontal direction')

    % Visualize line profile
    figure;
    plot(line_profile, 'LineWidth', 1.5);
    grid on; axis square;
    title(['Line profile through centre (' num2str(diameters(i)) ' mm)']);
    xlabel('Pixel index');
    ylabel('Pixel value - attenuation \mu');
end 

%% 3. AVG AND STD OF PLEXI ATTENUATION COEFFICIENT (plexi-only sections)

% Define a small square Region Of Interest around the centre (e.g. 20x20 pixels)
half_width = 10;  % ROI will be (2*half_width + 1) x (2*half_width + 1)

for i = 1:3   %indexes of 15 mm, 35 mm, 55 mm plexi-only sections

    avg_image_plexi = all_avg_images{i};

    cy = centers(i,1);
    cx = centers(i,2);

    row_min = max( round(cy - half_width), 1 );
    row_max = min( round(cy + half_width), dim_y );
    col_min = max( round(cx - half_width), 1 );
    col_max = min( round(cx + half_width), dim_x );

    % Extract properly centered ROI:
    roi_mu = avg_image_plexi(row_min:row_max, col_min:col_max);

    % Compute mean and standard deviation of μ in the ROI
    mu_mean = mean(roi_mu(:));

    % Display results
    fprintf('--- %d mm pure-plexi section ---\n', diameters(i));
    fprintf('Mean attenuation μ = %.4e [mm^-1]\n\n', mu_mean);

    figure;
    imagesc(avg_image_plexi); axis image; colormap gray; colorbar;
    hold on;
    rectangle('Position', [row_min, col_min, col_max-col_min, row_max-row_min], ...
        'EdgeColor', 'r', 'LineWidth', 1.5);
    plot(cy, cx, 'b+', 'MarkerSize', 15, 'LineWidth', 2);
    title([num2str(diameters(i)) ' mm plexi section with central ROI']);
    xlabel('Y - vertical direction'); ylabel('X - horizontal direction');

end 

% when the loop is over, the last values of roi_mu and cy, cx refer to 55
% mm sectionclear 
mu_std  = std(roi_mu(:));
SNR_plexi = mu_mean / mu_std;

fprintf('Std(μ)  = %.4e [mm^-1]\n', mu_std);
fprintf('SNR     = %.2f\n\n', SNR_plexi);


%% 4. LINE PROFILE THROUGH SECTION WITH 2 RODS

mm55_2rods_img = all_avg_images{4};  % section with two rods

% Centers of teflon and air rods
%xair = 273;     yair = 167;
%xteflon = 129;  yteflon = 208;

% Creates masks with air and teflon rods only to get their centers with
% circleCenterMask function
slices_to_average_uint16 = img_uint16(:, :, slices(4):(slices(4) + n_slices));
avg_image_uint16 = mean(slices_to_average_uint16,3);

air_mask = (avg_image_uint16 < 15000) & (avg_image_uint16 > 10000);
[x_air, y_air] = circleCenterMask(air_mask);

teflon_mask = (avg_image_uint16 > 27000);
[x_teflon, y_teflon] = circleCenterMask(teflon_mask);

%%

% Line equation model
m = (y_teflon - y_air) / (x_teflon - x_air);
q = y_air - m * x_air;
x = 1:dim_x;
y = m * x + q;
y_round = round(y);

% Valid values inside the image
valid = y_round>=1 & y_round<=dim_x & x>=1 & x<=dim_y;
yv = y_round(valid); 
xv = x(valid);

figure;
imagesc(mm55_2rods_img); axis image; colormap gray; colorbar;
hold on;
plot(xv, yv, 'r-', 'LineWidth', 2);
plot(x_air, y_air, 'gx', 'MarkerSize', 8, 'LineWidth', 2);
plot(x_teflon, y_teflon, 'bx', 'MarkerSize', 8, 'LineWidth', 2);
title('Line through teflon and air rods');
xlabel('Y - vertical direction'); ylabel('X - horizontal direction');

% Initialize full-length profile (one value per x); NaN where line is out of image
atten_x = nan(1, dim_x);

% Sample attenuation on the line for each valid x
idx = sub2ind(size(mm55_2rods_img), yv, xv);     % (row, col)
atten_x(xv) = mm55_2rods_img(idx);

% Line profile
figure;
plot(1:dim_x, atten_x, 'LineWidth', 1.5);
grid on;
xlabel('Pixel index');
ylabel('Pixel value - attenuation \mu');
title('Profile along line crossing air and teflon rods');

%% 5. ATTENUATION COEFFICIENTS OF PLEXI AND TEFLON IN SECTION WITH TWO RODS

col_min_tef = max( round(x_teflon - half_width), 1 );
col_max_tef = min( round(x_teflon + half_width), dim_y );
row_min_tef = max( round(y_teflon - half_width), 1 );
row_max_tef = min( round(y_teflon + half_width), dim_x );

% Extract properly centered ROI for plexi and teflon
roi_mu_plexi = mm55_2rods_img(row_min:row_max, col_min:col_max);
roi_mu_tef = mm55_2rods_img(row_min_tef:row_max_tef, col_min_tef:col_max_tef);

% Compute mean and standard deviation of μ in the ROI
mu_mean_plexi = mean(roi_mu_plexi(:));
mu_mean_tef = mean(roi_mu_tef(:));

% Display results for attenuation coefficient
fprintf('--- 55 mm plexi section with 2 rods---\n');
fprintf('Mean attenuation μ (plexi) = %.4e [mm^-1]\n\n', mu_mean_plexi);

fprintf('--- 55 mm plexi section with 2 rods---\n');
fprintf('Mean attenuation μ (teflon) = %.4e [mm^-1]\n', mu_mean_tef);


% Visualize ROI
figure;
imagesc(mm55_2rods_img); axis image; colormap gray; colorbar; hold on;

rectangle('Position', [row_min, col_min, col_max-col_min, row_max-row_min], ...
          'EdgeColor', 'r', 'LineWidth', 1.5);
rectangle('Position', [col_min_tef, row_min_tef, row_max_tef-row_min_tef, col_max_tef-col_min_tef], ...
          'EdgeColor', 'r', 'LineWidth', 1.5);

plot(x_teflon, y_teflon, 'b+', 'MarkerSize', 15, 'LineWidth', 2);
plot(cy, cx, 'b+', 'MarkerSize', 15, 'LineWidth', 2);

title('55 mm plexi section with two rods (ROIs)');
xlabel('Y - vertical direction'); ylabel('X - horizontal direction');

diary off

%% Save figures

outDir = uigetdir(pwd, 'Select folder to save images');
figs = findall(0,'Type','figure');

for k = 1:length(figs)
    exportgraphics(figs(k), fullfile(outDir, sprintf('Figure_%02d.png', k)));
end