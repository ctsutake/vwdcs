-----------------------------------------------------------------------------

        Vaguelette-wavelet deconvolution via compressive sampling            
 
-----------------------------------------------------------------------------

Written by  : Chihiro Tsutake
Affiliation : University of Fukui
E-mail      : ctsutake@icloud.com
Created     : April 2019

-----------------------------------------------------------------------------
    Contents
-----------------------------------------------------------------------------

demo.m      : Main algorithm file
img.m       : Read ideal image 
psf.m       : Read blur impulse response
breg.m      : Split-Bregman method for our CS-based recovery
cyclespin.m : Cycle-spinning for our modified SBT-based denoising
denoise.m   : Modified SBT-based denoising
img/        : Original PNG images
mex/        : Mex sources (https://statweb.stanford.edu/~wavelab/)

-----------------------------------------------------------------------------
    Preliminary
-----------------------------------------------------------------------------

1) `cd mex'
2) `mex FWT2_PO.c' and `mex IWT2_PO.c'

-----------------------------------------------------------------------------
    Usage
-----------------------------------------------------------------------------

1) Change parameters in demo.m and ARGs of the functions img() and psf().
2) Running `demo' generates the following images.

    -- obs.png (observed image)
    -- den.png (denoised image)
    -- res.png (restored image)

-----------------------------------------------------------------------------
    Feedback
-----------------------------------------------------------------------------

If you have any questions, please contact me.
