% cyclespin.m
%
% Cycle-spinning
%
% x - input image
% i - vertical shift
% j - horizontal shift
%
% Written by  : Chihiro Tsutake
% Affiliation : University of Fukui
% E-mail      : ctsutake@icloud.com
% Created     : April 2019
%

function s = cyclespin(x, i, j)

    [l1, l2] = size(x);
    z = x((l1 + 1 - i):l1, :);
    z(i + 1:l1, :) = x(1:(l1 - i), :);
    s = z(:, (l2 + 1 - j):l2);
    s(:, j + 1:l2) = z(:, 1:(l2 - j));
