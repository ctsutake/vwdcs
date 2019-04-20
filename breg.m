% breg.m
%
% Split-bregman method for Eq.(10)
%
% OMEGA - Fourier-domain indices in Eq.(8)
% X_hat - partial DFT coefficients in Eq.(9)
% mu    - regularization parameter in Eq.(10)
% eps   - convergence criterion in Eq.(11)
%
% Written by  : Chihiro Tsutake
% Affiliation : University of Fukui
% E-mail      : ctsutake@icloud.com
% Created     : April 2019
%

function x_hat = breg(OMEGA, X_hat, mu, eps)

    x_hat = zeros(256);

    % Variables
    prm_u = randn(256);
    aux_s = zeros(256);
    aux_t = zeros(256);
    lag_v = zeros(256);
    lag_w = zeros(256);

    % Kernel (Fourier domain)
    ker = zeros(256);
    ker([end, 1, 2], [end, 1, 2]) = [0, -1, 0; -1, 4, -1; 0, -1, 0];
    ker = real(fft2(ker));
    ker(OMEGA) = ker(OMEGA) + 1;

    for rho = 1:65536
        % Copy old
        old_u = prm_u;

        % Update primal variable
        rhs = aux_s - lag_v; prm_u = -rhs + [rhs(:, 256), rhs(:, 1:255)];
        rhs = aux_t - lag_w; prm_u = -rhs + [rhs(256, :); rhs(1:255, :)] + prm_u;
        prm_u = real(ifft2((X_hat + fft2(prm_u)) ./ ker));

        % Update auxiliary variable
        tmp1 = -prm_u + [prm_u(:, 2:256), prm_u(:, 1)];
        tmp2 = -prm_u + [prm_u(2:256, :); prm_u(1, :)];
        tmp3 = sqrt((tmp1 + lag_v).^2 + (tmp2 + lag_w).^2);
        aux_s = max(tmp3 - mu, 0) .* (tmp1 + lag_v) ./ tmp3;
        aux_t = max(tmp3 - mu, 0) .* (tmp2 + lag_w) ./ tmp3;

        % Update Lagrangian multiplier
        lag_v = lag_v + tmp1 - aux_s;
        lag_w = lag_w + tmp2 - aux_t;

        % Convergence criterion
        if norm(prm_u(:) - old_u(:)) / norm(old_u(:)) < eps
            break;
        end

    end

    x_hat = prm_u;
