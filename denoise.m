% denoise.m
%
% Modified SBT-based denoising in Sect. III-A
%
% noisy - noisy image
% sgm   - standard deviation of noise
% dlt   - maximum shift value
%
% Written by  : Chihiro Tsutake
% Affiliation : University of Fukui
% E-mail      : ctsutake@icloud.com
% Created     : April 2019
%

function out = denoise(noisy, sgm, dlt)

    addpath('mex');
    out = zeros(256);

    % Symmlet6 reconstruction filter
    qmf = [-0.007800708325034148, 0.0017677118642428036, 0.04472490177066578, -0.021060292512300564, -0.07263752278646252, 0.3379294217276218, 0.787641141030194, 0.4910559419267466, -0.048311742585633, -0.11799011114819057, 0.0034907120842174702, 0.015404109327027373];

    % Constants in Eq.(2)
    C = 4.50524 * sgm^2;
    L = 3;
    j0 = 1;

    % Kernels
    k1 = ones(L) / L^2;
    k2 = ones(L);

    for dn1 = 0:dlt

        for dn2 = 0:dlt
            % Forward cycle-spinning
            noistrans = cyclespin(noisy, dn1, dn2);

            % Forward wavelet transform
            wc = FWT2_PO(noistrans, j0, qmf);

            % Stein block thresholding (block average is computed via filtering and down/up sampling)
            % j = 7, i = 1, 2, 3
            thr = wc(001:128, 129:256).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(001:128, 129:256) = max(1 - C * thr.^(-1), 0) .* wc(001:128, 129:256);
            thr = wc(129:256, 001:128).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(129:256, 001:128) = max(1 - C * thr.^(-1), 0) .* wc(129:256, 001:128);
            thr = wc(129:256, 129:256).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(129:256, 129:256) = max(1 - C * thr.^(-1), 0) .* wc(129:256, 129:256);
            % j = 6, i = 1, 2, 3
            thr = wc(001:064, 065:128).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(001:064, 065:128) = max(1 - C * thr.^(-1), 0) .* wc(001:064, 065:128);
            thr = wc(065:128, 001:064).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(065:128, 001:064) = max(1 - C * thr.^(-1), 0) .* wc(065:128, 001:064);
            thr = wc(065:128, 065:128).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(065:128, 065:128) = max(1 - C * thr.^(-1), 0) .* wc(065:128, 065:128);
            % j = 5, i = 1, 2, 3
            thr = wc(001:032, 033:064).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(001:032, 033:064) = max(1 - C * thr.^(-1), 0) .* wc(001:032, 033:064);
            thr = wc(033:064, 001:032).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(033:064, 001:032) = max(1 - C * thr.^(-1), 0) .* wc(033:064, 001:032);
            thr = wc(033:064, 033:064).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(033:064, 033:064) = max(1 - C * thr.^(-1), 0) .* wc(033:064, 033:064);
            % j = 4, i = 1, 2, 3
            thr = wc(001:016, 017:032).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(001:016, 017:032) = max(1 - C * thr.^(-1), 0) .* wc(001:016, 017:032);
            thr = wc(017:032, 001:016).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(017:032, 001:016) = max(1 - C * thr.^(-1), 0) .* wc(017:032, 001:016);
            thr = wc(017:032, 017:032).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(017:032, 017:032) = max(1 - C * thr.^(-1), 0) .* wc(017:032, 017:032);
            % j = 3, i = 1, 2, 3
            thr = wc(001:008, 009:016).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(001:008, 009:016) = max(1 - C * thr.^(-1), 0) .* wc(001:008, 009:016);
            thr = wc(009:016, 001:008).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(009:016, 001:008) = max(1 - C * thr.^(-1), 0) .* wc(009:016, 001:008);
            thr = wc(009:016, 009:016).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 2, 1:end - 2) = kron(avg, k2); wc(009:016, 009:016) = max(1 - C * thr.^(-1), 0) .* wc(009:016, 009:016);
            % j = 2, i = 1, 2, 3
            thr = wc(001:004, 005:008).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(001:004, 005:008) = max(1 - C * thr.^(-1), 0) .* wc(001:004, 005:008);
            thr = wc(005:008, 001:004).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(005:008, 001:004) = max(1 - C * thr.^(-1), 0) .* wc(005:008, 001:004);
            thr = wc(005:008, 005:008).^2; avg = conv2(thr, k1, 'valid'); avg = avg(1:L:end, 1:L:end); thr(1:end - 1, 1:end - 1) = kron(avg, k2); wc(005:008, 005:008) = max(1 - C * thr.^(-1), 0) .* wc(005:008, 005:008);
            % j = 1, i = 1, 2, 3
            wc(001:002, 003:004) = max(1 - C * wc(001:002, 003:004).^(-2), 0) .* wc(001:002, 003:004);
            wc(003:004, 001:002) = max(1 - C * wc(003:004, 001:002).^(-2), 0) .* wc(003:004, 001:002);
            wc(003:004, 003:004) = max(1 - C * wc(003:004, 003:004).^(-2), 0) .* wc(003:004, 003:004);

            % Inverse wavelet transform
            dout = IWT2_PO(wc, j0, qmf);

            % Inverse cycle-spinning
            dout = cyclespin(dout, 256 - dn1, 256 - dn2);

            out = out + dout;
        end

    end

    % Averaging
    out = out / ((dlt + 1) * (dlt + 1));
