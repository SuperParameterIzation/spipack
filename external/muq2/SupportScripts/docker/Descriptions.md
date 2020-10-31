## Descriptions

### `muq-build`
Defines a docker image based on `debian:latest` with all of MUQ's dependencies installed.  This container does not have MUQ installed, but is used in our continuous integration pipeline for testing.

### `muq-hippylib`
Defines a docker image with installations of MUQ, Fenics, and [hIPPYlib](https://hippylib.github.io/).  This is based on the `quay.io/fenicsproject/stable` docker image.
