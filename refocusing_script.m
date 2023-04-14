clear;


%% parameters

% parameters of microscope setup
NA = .4;
lambda = 633E-9;
k = 2 * pi / lambda;

% image stretch parameters
symm_image = 1;             % pixel aspect ratio y/x
yoffs = -10;                % apply shift in y when cropping image back to square aspect ratio

% refocusing parameters;
cutoff_k = k*NA*1.5;        % upper limit for k, empirically chosen

% line out parameters
flag_amplitude = 0;         % specify whether final image is treated as intensity or amplitude
crop_xy = [21,680,21,680];  % crop image after refocusing



%% get image
% fname = 'reco_z0'; load(fname); delta_z = 0*1e-6; symm_image  = 60/83; flag_amplitude = 1;
% fname = 'reco_z20'; load(fname);  delta_z = -16.0*1e-6; symm_image = 60/83; flag_amplitude = 1; flag_amplitude = 1;
fname = 'reco_z40';  load(fname);  delta_z = -40*1e-6; symm_image  = 62/83; linex = 157; linex2 = 540; flag_amplitude  = 1;
% fname = 'reco_z60'; load(fname);  delta_z = -61*1e-6; symm_image  = 63.5/83; flag_amplitude = 1;
% fname = 'reco_zm20'; load(fname); delta_z = 23.5*1e-6; symm_image  = 64/83;  flag_amplitude = 1;

% fname = 'noninter_z40'; load(fname); jpk_reco = imrotate(jpk_reco, 1.3,'bilinear'); jpk_reco=jpk_reco(1:720,35:35+719); delta_z = 0; symm_image  = 62.4/83; yoffs = -12;
% fname = 'noninter_z0'; load(fname); jpk_reco = imrotate(jpk_reco, 1.3,'bilinear'); jpk_reco=jpk_reco(1:720,35:35+719);  delta_z = 0; symm_image  = 62.4/83; yoffs = -12;  linex = 170; linex2 = 553;


%% pre-process image

% rotate image
mat_image = rot90(jpk_reco); % imrotate(jpk_reco,90);
[Ny Nx] = size(mat_image);

% stretch image in y to achieve pixel aspect ratio of 1 (square pixel), then crop to Nx
if(symm_image<1)
    mat_image = imresize(mat_image, [round(Ny/symm_image) Nx]);
    [Ny Nx] = size(mat_image);

    mat_image = mat_image( round(Ny/2-Nx/2)+yoffs:round(Ny/2-Nx/2)+Nx-1+yoffs , :);
    [Ny Nx] = size(mat_image);
end


% apply window function to image to pull signal to mean value towards the edges (ensures proper FFT)
m = mean(mean(mat_image));
window1d = tukeywin(Ny,.1); mat_window = window1d .* window1d';
image_oof = mat_image.*mat_window + m*(1-mat_window);


%% preparations for refocusing

% define sampling rate
fs = Ny/80E-6*0.98;
dx = 1/fs;
dy = 1/fs;

% calculate matrix for kz(k_||/2)
kx_max = 2*pi/(2*dx);
ky_max = 2*pi/(2*dy);
kx = linspace(-kx_max,kx_max,Nx);
ky = linspace(-ky_max,ky_max,Ny);
[KY, KX] = meshgrid(ky,kx);
k_pll = sqrt(KX.^2+KY.^2);
k_mask = (abs(KX)<k) .* (abs(KY)<k);
kz = sqrt(k^2-(k_pll/2).^2);

% calculating phaseing matrix
refocusing_op = 1./(kz/k) .* exp( 2i * delta_z * kz );
refocusing_op( k_pll > k ) = 0; % remove evanescent components


%% perform refocusing
fft_data = fftshift(fft2( image_oof ));
fft_data_refocused = fft_data .* refocusing_op;
mat_refocused = ifft2(fftshift(fft_data_refocused));


%% post process

% set phase on substrate to zero
section_image = mat_refocused(:, [1:50 668:Nx]);
avg_phase = mean( mean( section_image./abs(section_image), 2 ), 1);
mat_refocused = mat_refocused * conj(avg_phase);



%% plot
figure(1); clf;
subplot(3,1,1); imagesc(abs(mat_image)); title('Original image, abs');
subplot(3,1,2); imagesc(angle(mat_image)); title('Original image, angle');
subplot(3,1,3); imagesc(abs(mat_window)); title('Window function');
colormap 'gray';

figure(2); clf;
subplot(3,1,1); imagesc( kx,ky,log(abs(fft_data)) ); title('FT of image, abs');
subplot(3,1,2); imagesc( kx,ky,abs(refocusing_op)); title('Refocusing operator, abs');
subplot(3,1,3); imagesc( kx,ky,angle(refocusing_op)); title('Refocusing operator, angle');
colormap 'gray';

figure(3); clf;
subplot(2,1,1); imagesc( abs(mat_refocused) ); title('Refocused image, abs');
subplot(2,1,2); imagesc( angle(mat_refocused) ); title('Refocused image, angle');
colormap 'gray';

%% plot line outs

% crop image to remove empty space at border
if(exist('crop_xy','var'))
    mat_image = mat_image(crop_xy(1):crop_xy(2), crop_xy(3):crop_xy(4));
    mat_refocused = mat_refocused(crop_xy(1):crop_xy(2), crop_xy(3):crop_xy(4));
end
[Ny Nx] = size(mat_image);


if exist('linex','var')

    % take vertical line outs
    lineoutx = mat_refocused(:,linex);
    lineoutx2 = mat_refocused(:,linex2);
    if flag_amplitude
        % this is for the SOH data which is in amplitude and needs to be
        % squared for direct comparison with intensity
        lineoutx = abs(lineoutx).^2;
        lineoutx2 = abs(lineoutx2).^2;
    end

    figure(4); clf;
    subplot(2,1,1);plot(abs(lineoutx)); title('Line out group 8');
    subplot(2,1,2);plot(abs(lineoutx2));  title('Line out group 9');
    xarray = linspace(0, 46.874/(569-152)*(Ny-1),Ny); % 46.874um = distance between left most and right most bar, 569-152 distance in pixel
%     dlmwrite( strcat(fname,'_dz',num2str(delta_z*1e6),'_linex2_',num2str(linex2),'.txt'), [transpose(xarray) abs(lineoutx2).^2 angle(lineoutx2)], ' ');

end
