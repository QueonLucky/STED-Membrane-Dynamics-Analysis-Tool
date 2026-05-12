clear; close all; clc;

filename = 'membrane.tif';  
time_interval = 0.5;                 
pixel_size = 0.06;                   % 60nm

%% 

fprintf('readingstack...\n');
info = imfinfo(filename);
T = numel(info);                     
Y = info(1).Height;                 
X = info(1).Width;                   

I = zeros(Y, X, T, 'single');

for t = 1:T
    I(:,:,t) = im2single(imread(filename, t));
    if mod(t, 50) == 0
        fprintf('alreadyread %d / %d 帧\n', t, T);
    end
end
fprintf('Donereadfor %d x %d x %d\n', Y, X, T);


fprintf('Estimateandsubtractthebackground...\n');

sorted_I = sort(I(:), 'ascend');
bg_cut = round(0.05 * numel(sorted_I));
bg_value = mean(sorted_I(1:bg_cut));
fprintf('thebackground：%.4f\n', bg_value);


I_bgsub = max(I - bg_value, 0);
meanIntensity = squeeze(mean(I_bgsub, [1 2]));  
time_vec = (0:T-1)' * time_interval;            

fprintf('Fittingofphotobleachingtrends...\n');
f = fit(time_vec, meanIntensity, 'exp2');
bleach_trend = f(time_vec);   
figure;
plot(time_vec, meanIntensity, 'b-', time_vec, bleach_trend, 'r--');
xlabel('Time (s)'); ylabel('MeanIntensity'); legend('RawData','photobleachingtrends');
title('Fittingofphotobleachingtrends');

I_corrected = I_bgsub ./ reshape(bleach_trend, 1, 1, T);
I_corrected(~isfinite(I_corrected)) = 0;

fprintf('Apply3DGaussianfilteringfornoisereduction...\n');
I_smoothed = imgaussfilt3(I_corrected, [0.8, 0.8, 0.5]);

%% 

global_signal = squeeze(mean(I_smoothed, [1, 2]));  
global_signal_detrend = detrend(global_signal);    

figure;
subplot(2,1,1);
plot(time_vec, global_signal); xlabel('Time (s)'); ylabel('Intensity');
title('Global Average Fluorescence Intensity');
subplot(2,1,2);
plot(time_vec, global_signal_detrend); xlabel('Time (s)'); ylabel('Intensity');
title('Fluctuations after detrending');

%% 

Fs = 1 / time_interval;  
L = length(global_signal_detrend);
[Pxx, f] = periodogram(global_signal_detrend, [], L, Fs);

figure;
plot(f, Pxx);
xlabel('Frequency (Hz)'); ylabel('Power Spectral Density');
xlim([0 Fs/2]);   
title('Global Signal Power Spectrum');

[~, idx] = max(Pxx(f > 0.01)); 
peak_freq = f(idx);
fprintf('The main oscillation frequency is：%.4f Hz (Frequency %.2f s)\n', peak_freq, 1/peak_freq);


%% 

pixel_std = std(I_smoothed, 0, 3);  

figure;
imagesc(pixel_std);
axis image; colorbar;
title('Fluorescence Intensity Fluctuation Amplitude(Standard Deviation)');
xlabel('x (pixel)'); ylabel('y (pixel)');

%% 

pixel_mean = mean(I_smoothed, 3);
pixel_cv = pixel_std ./ (pixel_mean + 1e-6);  

figure;
imagesc(pixel_cv);
axis image; colorbar;
title('Fluorescence Intensity Fluctuation Amplitude(CV)');

%% 

roi_radius = 10;   
center1 = [50, 50];
center2 = [50, 100];  


y1 = (center1(1)-roi_radius):(center1(1)+roi_radius);
x1 = (center1(2)-roi_radius):(center1(2)+roi_radius);
signal1 = squeeze(mean(I_smoothed(y1, x1, :), [1, 2]));

y2 = (center2(1)-roi_radius):(center2(1)+roi_radius);
x2 = (center2(2)-roi_radius):(center2(2)+roi_radius);
signal2 = squeeze(mean(I_smoothed(y2, x2, :), [1, 2]));


s1_det = detrend(signal1);
s2_det = detrend(signal2);


[corr, lags] = xcorr(s1_det, s2_det, 'normalized');
[~, max_idx] = max(corr);
time_lag = lags(max_idx) * time_interval;   
distance_px = abs(center2(2) - center1(2)); 
velocity_px_per_sec = distance_px / time_lag;
velocity_um_per_sec = velocity_px_per_sec * pixel_size;

fprintf('The wave speed is approximately %.2f pixel/s，as %.2f micrometers/s\n', ...
    velocity_px_per_sec, velocity_um_per_sec);

%% 

figure;
imagesc(I_smoothed(:,:,1));  
axis image; colormap gray;
title('Drawaline and doubleclick');

% Define ROI
h = drawfreehand('Color','red','LineWidth',2); 

position = h.Position;  
disp('your ROI：');
disp(position);


num_points = 100;  
T = size(I_smoothed, 3);
kymo_line = zeros(num_points, T, 'single');

fprintf('Drawing kymograph...\n');
for t = 1:T
    [cx, cy, c] = improfile(I_smoothed(:,:,t), position(:,1), position(:,2), num_points);
    kymo_line(:, t) = c;
end

figure;
imagesc(1:num_points, (0:T-1)*time_interval, kymo_line');
xlabel('Drawline'); ylabel('Time (s)');
title('KymographofROI');
colormap gray; colorbar;


kymo_for_radon = kymo_line';  
theta = -89:1:89;  
[R, xp] = radon(kymo_for_radon, theta);
[~, max_idx] = max(var(R, 0, 1));   
best_theta = theta(max_idx);
figure;
imagesc(theta, xp, R);
xlabel('angle (°)'); ylabel('Projection position');
title('Radon Transform'); colormap hot; colorbar;

%% 

spatial_res = 1;   
time_res = time_interval;   
if abs(best_theta) > 0.1
    velocity = tan(best_theta * pi/180) * (spatial_res / time_res);
    fprintf('Estimate wave velocity:%.2f pixel/s\n', velocity);
else
    fprintf('The stripes are nearly horizontal, with extremely low or stationary wave speed\n');
end

%% 

[Gx, Gt] = gradient(kymo_line);
angle_map = atan2(Gt, Gx);
histogram(angle_map(:), -pi/2:0.05:pi/2);
xlabel('angle (rad)'); ylabel('Frequency');
title('Local Stripe Orientation Distribution');
%% 


opticFlow = opticalFlowFarneback;
velocity_x = zeros(Y, X, T-1, 'single');
velocity_y = zeros(Y, X, T-1, 'single');
for t = 2:T
    frame_current = I_smoothed(:,:,t);
    frame_prev = I_smoothed(:,:,t-1);
    flow = estimateFlow(opticFlow, frame_prev);
    velocity_x(:,:,t-1) = flow.Vx;   
    velocity_y(:,:,t-1) = flow.Vy;
end

speed_map = sqrt(mean(velocity_x.^2 + velocity_y.^2, 3));
imagesc(speed_map); colorbar; title('Farneback');


%% 

I = single(I_smoothed);
[Y, X, T] = size(I);
I = I - mean(I, [1 2]);
I_detrended = zeros(size(I), 'single');
for y = 1:Y
    for x = 1:X
        I_detrended(y,x,:) = detrend(squeeze(I(y,x,:)), 'linear');
    end
end

win_x = hann(X)';      
win_y = hann(Y);       
win_t = reshape(hann(T), 1, 1, T);

win3d = win_y .* win_x .* win_t;
I_windowed = I_detrended .* win3d;

F = fftn(I_windowed);
F_shift = fftshift(F);
P = abs(F_shift).^2;
P_log = log10(P + 1);   

dx = pixel_size;      
dy = pixel_size;
dt = time_interval;   

kx_axis = (-floor(X/2):ceil(X/2)-1) / (X*dx);
ky_axis = (-floor(Y/2):ceil(Y/2)-1) / (Y*dy);
f_axis  = (-floor(T/2):ceil(T/2)-1) / (T*dt);

[~, ky0_idx] = min(abs(ky_axis));
slice_kx_f = squeeze(P_log(ky0_idx, :, :));   % [X, T]

figure;
center_val = max(P_log(:));
bg_level = mean(P_log(P_log < prctile(P_log(:), 90)));
clim_low = bg_level + 0.5 * (center_val - bg_level);   
clim_high = center_val;

imagesc(kx_axis, f_axis, slice_kx_f');
xlabel('k_x (1/\mum)'); ylabel('f (Hz)');
title('ω-k spectrum (k_y=0 section)');
colorbar; axis xy; clim([clim_low, clim_high]);  


%% 

frame_mid = round(T/2);
corr_map = normxcorr2(I_smoothed(:,:,frame_mid), I_smoothed(:,:,frame_mid));
imagesc(corr_map); axis image;

%% 
global_signal = squeeze(mean(I_smoothed, [1 2]));
[Pxx, f] = periodogram(detrend(global_signal), [], [], 1/time_interval);
plot(f, Pxx);
xlabel('Frenquency (Hz)'); ylabel('Power');
title('Spatially Averaged Signal Spectrum');


%%

data_matrix = reshape(I_smoothed, Y*X, T)';
data_centered = data_matrix - mean(data_matrix, 1);

[coeff, score, latent] = pca(data_centered);

figure;
bar(cumsum(latent)/sum(latent)*100);
xlabel('PCA'); ylabel('Cumulative contribution(%)');
title('PCA cumulative explained variance');

figure;
for k = 1:4
    mode_map = reshape(coeff(:,k), Y, X);
    subplot(2,2,k);
    imagesc(mode_map); axis image; colorbar;
    title(sprintf('PC %d (%.1f%%)', k, latent(k)/sum(latent)*100));
end

