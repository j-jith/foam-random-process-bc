# randomProcessFixedValue

This is an OpenFOAM boundary condition for generating random gusts at a
velocity inlet.

**This is a work in progress. Expect breaking changes.**

## Usage

The random gust is generated using the Spectral Representation Method from a
specified power spectral density (PSD). At the moment, this boundary condition
can generate a horizontal/vertical spatially uniform random gust.

The following parameters are used to specify the gust:

- `gustDir`: a vector representing the direction of the gust
- `scaleFac`: a scaling factor for the amplitude of the gust
- `circFreq`: a list of frequencies at which the PSD is given
- `psd`: a list representing the PSD at each frequency in `circFreq`

## Example

Two example cases are provided in the [run](run/) directory:

- [rectangularCavity](run/rectangularCavity/): A rectangular cavity with a
  velocity inlet
- [rectangularCavityPar](run/rectangularCavityPar): Same rectangular cavity
  modified to run in parallel

Please have a look at the [0/U](run/rectangularCavityPar/0/U) file to
understand how the gust is to be specified.
