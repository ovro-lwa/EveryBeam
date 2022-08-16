External dependencies are included via git submodules, see .gitmodules in the root of the source.

Please note that the dependencies are fixed to a specific commit, to make sure tests do not break as a 
result of silently updating the submodules. Dependencies are checked out on the following commits:

- `aocommon`: dd3722b9c3c3b8f0eb14f75db4b0bfb80cc91a68 (https://gitlab.com/aroffringa/aocommon)
- `eigen`: `3.3.7` (https://gitlab.com/libeigen/eigen/-/tags/3.3.7)
- `pybind11`: `v2.5.0` (https://github.com/pybind/pybind11/releases/tag/v2.5.0)

The file npy.hpp comes from the libnpy library that can be found here:

- https://github.com/llohse/libnpy

See that file for its license.