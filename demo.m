% demo.m (256x256 grayscale image is only supported)
%
% Vaguelette-wavelet deconvolution via compressive sampling
%
% Written by  : Chihiro Tsutake
% Affiliation : University of Fukui
% E-mail      : ctsutake@icloud.com
% Created     : April 2019
%

% Standard deviation of Gaussian noise in Eq.(1)
sgm = 1;
% Maximum shift value in Eq.(7)
dlt = 3;
% Threshold in Eq.(8)
tau = 0.08;
% Regularization parameter in Eq.(10)
mu = 1;
% Convergence criterion in Eq.(11)
eps = 1E-4;

% The functions img(1), img(2), ..., img(6) read Pillars, Flowers, Swans_1, Sphynx, Bush, and Friends_3 in Fig.4, reps..
x_idl = img(3);

% The functions psf(1), psf(2), ..., psf(6) read h_a, h_b, ..., h_f in Fig.5, reps..
h_pad = psf(5);

% Observation based on Eq.(1)
H_pad = fft2(h_pad);
X_idl = fft2(x_idl);
Y_idl = X_idl .* H_pad;
y_idl = real(ifft2(Y_idl));
y_obs = y_idl + randn(256) * sgm;

% Modified SBT-based denoising in Sect. III-A
y_hat = denoise(y_obs, sgm, dlt);

% Set OMEGA in Eq.(8)
OMEGA = find(abs(H_pad) >= tau);

% Partial inversion in Eq.(9)
X_hat = zeros(256); X_hat(OMEGA) = fft2(y_hat)(OMEGA) ./ H_pad(OMEGA);

% Split-Bregman method for Eq.(10)
x_hat = breg(OMEGA, X_hat, mu, eps);

% Write images
imwrite(uint8(y_obs), 'obs.png'); % observed
imwrite(uint8(y_hat), 'den.png'); % denoised
imwrite(uint8(x_hat), 'res.png'); % restored

% PSNR
printf("PSNR(x_idl, x_hat) = %4.2f\n", psnr(uint8(x_idl), uint8(x_hat)));
