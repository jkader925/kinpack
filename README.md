# kinpack
Stellar kinematic analysis for observations and simulations.

This tool was written in Interactive Data Language (IDL) to reconstruct 2-D resolved velocity and velocity dispersion fields of galaxies. The original dataset consists of radial velocities of resolved globular clusters (GCs) in nearby giant quenched galaxies, obtained with the DEIMOS fiber-fed multi-object spectrograph on the Keck II telescope on Mauna Kea. Kinematic map construction is based on the geostatistics technique known as "Kriging", essentially interpolating a smooth value map based on measured values at discrete positions. Kriging is based on Gaussian process regression, and is also known as Wiener-Kolmogorov prediction.

See https://www.justinkader.com/research "Stellar Kinematics in Early-type Galaxies" tab for example output.
