% img.m
%
% Read ideal image utilized in our paper
%
% Written by  : Chihiro Tsutake
% Affiliation : University of Fukui
% E-mail      : ctsutake@icloud.com
% Created     : April 2019
%

function x = img(t)

    switch t
        case 1
            x = double(imread('./img/pillars.png'));
        case 2
            x = double(imread('./img/flowers.png'));
        case 3
            x = double(imread('./img/swans_1.png'));
        case 4
            x = double(imread('./img/sphynx.png'));
        case 5
            x = double(imread('./img/bush.png'));
        case 6
            x = double(imread('./img/friends_3.png'));
    end
