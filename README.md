# STED-Membrane-Dynamics-Analysis-Tool
A MATLAB-based image processing pipeline for analyzing fluorescence dynamics in STED microscopy videos of cell membranes. The tool extracts temporal fluctuations, wave propagation speeds, and spatial activity patterns from time-lapse stacks.

Features
Automated preprocessing – background subtraction, photobleaching correction, and 3D Gaussian filtering

Global dynamics analysis – average intensity trends, detrended fluctuations, and power spectrum (periodogram)

Spatial activity mapping – standard deviation and coefficient of variation (CV) maps of intensity fluctuations

Wave speed estimation – via cross-correlation between ROIs, Radon transform of kymographs, and local gradient orientation

Optical flow – Farneback method for pixel-wise velocity fields

3D Fourier analysis – frequency–wavenumber (ω-k) spectrum to reveal dispersion relations

PCA decomposition – identify dominant spatiotemporal modes with explained variance

Requirements
MATLAB R2019b or later

Image Processing Toolbox

Signal Processing Toolbox

Statistics and Machine Learning Toolbox

Computer Vision Toolbox (for optical flow)

Input
A multi-frame TIFF file (e.g., membrane.tif) containing a time-lapse STED image stack.
Parameters to set inside the script:

time_interval – time between frames (seconds)

pixel_size – physical pixel size (µm)

Output
Figures generated during execution:

Photobleaching fit and corrected intensity curve

Global fluorescence intensity and detrended fluctuations

Power spectral density (dominant frequency)

Spatial maps of fluctuation amplitude (std and CV)

Cross-correlation lag and derived wave speed

Kymograph along a user-drawn line

Radon transform for kymograph angle analysis

Local stripe orientation histogram

Farneback optical flow speed map

ω-k spectrum (kx vs frequency slice)

PCA cumulative variance and top spatial modes

Command-line outputs include background intensity, dominant oscillation frequency, and estimated wave propagation speed (pixel/s and µm/s).

Usage
Place your STED time-lapse TIFF in the same folder as the script.

Adjust filename, time_interval, and pixel_size at the top of membrane_analysis_for_STED.m.

Run the script in MATLAB.

When prompted, draw a freehand line on the displayed image and double-click to generate a kymograph.

Example
matlab
filename = 'membrane.tif';
time_interval = 0.5;   % 2 fps
pixel_size = 0.06;     % 60 nm
Notes
The script assumes the membrane signal is brighter than background.

Photobleaching is fitted with a double-exponential (exp2) model.

For the Radon-based velocity estimation, the kymograph is assumed to have diagonal stripes whose angle corresponds to wave speed.

License
MIT – feel free to use and modify for your research.

Author
Developed for STED membrane dynamics analysis.
