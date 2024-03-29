% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rcpp_simulateFullWaveform.R
\name{simulateFullWaveform}
\alias{simulateFullWaveform}
\title{Simulate full wave form from point cloud LAS file}
\usage{
simulateFullWaveform(input, output, coords, inList = FALSE,
  ground = FALSE, hdf = TRUE, waveID = character(0), pSigma = -1,
  pFWHM = 15, readPulse = character(0), fSigma = 5.5,
  wavefront = character(0), res = 0.15, LVIS = FALSE,
  topHat = FALSE, sideLobe = FALSE, lobeAng = 0,
  checkCover = FALSE, maxScanAng = 1e+06, decimate = 1,
  pBuff = 0.2, maxBins = 1024L, countOnly = FALSE,
  pulseAfter = FALSE, pulseBefore = FALSE, noNorm = FALSE,
  noOctree = FALSE, octLevels = 0L, nOctPix = 40L, decon = FALSE,
  indDecon = FALSE, readWave = FALSE, listFiles = FALSE,
  keepOld = FALSE, useShadow = FALSE, polyGround = FALSE,
  nnGround = FALSE, seed = numeric(0))
}
\arguments{
\item{input}{name. lasfile input filename}

\item{output}{name. output filename}

\item{coords}{lon lat matrix. footprint coordinate in same system as lasfile}

\item{inList}{list. input file list (ASCII file) for multiple files}

\item{ground}{record separate ground and canopy waveforms}

\item{hdf}{write output as HDF5, default TRUE. Best with gridded or list of coords}

\item{waveID}{id.}

\item{pSigma}{sig. set Gaussian pulse width as 1 sigma}

\item{pFWHM}{fhwm. set Gaussian pulse width as FWHM in ns}

\item{readPulse}{file. read pulse shape and width from a file instead of making Gaussian}

\item{fSigma}{sig. set footprint width}

\item{wavefront}{file. read wavefront shape from file instead of setting Gaussian. Note that footprint width is still set by fSigma}

\item{res}{res. range resolution of waveform digitisation to output, in units of ALS data}

\item{LVIS}{use LVIS pulse length, sigma=6.25m}

\item{topHat}{use a top hat wavefront}

\item{sideLobe}{use side lobes}

\item{lobeAng}{ang. lobe axis azimuth}

\item{checkCover}{check that the footprint is covered by ALS data. Do not output if not}

\item{maxScanAng}{ang. maximum scan angle, degrees}

\item{decimate}{x. probability of accepting an ALS beam}

\item{pBuff}{s. point reading buffer size in Gbytes}

\item{maxBins}{Optional: for HDF5, limit number of bins to save trimming.}

\item{countOnly}{only use count method}

\item{pulseAfter}{apply the pulse smoothing after binning for computational speed, at the risk of aliasing (default)}

\item{pulseBefore}{apply the pulse smoothing before binning to avoid the risk of aliasing, at the expense of computational speed}

\item{noNorm}{don't normalise for ALS density}

\item{noOctree}{do not use an octree}

\item{octLevels}{n. number of octree levels to use}

\item{nOctPix}{n. number of octree pixels along a side for the top level}

\item{decon}{deconvolve}

\item{indDecon}{deconvolve individual beams}

\item{readWave}{read full-waveform where available}

\item{listFiles}{list files. Do not read them}

\item{keepOld}{do not overwrite old files, if they exist}

\item{useShadow}{account for shadowing in discrete return data through voxelisation}

\item{polyGround}{find mean ground elevation and slope through fitting a polynomial}

\item{nnGround}{find mean ground elevation and slope through nearest neighbour}

\item{seed}{n. seed number for random numbers generator}
}
\description{
Simulate full wave form from point cloud LAS file
}
