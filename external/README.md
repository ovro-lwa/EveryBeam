External dependencies are included via git submodules, see .gitmodules in the root of the source.

Please note that the dependencies are fixed to a specific commit, to make sure tests do not break due to
a silent updateing of the submodules. Dependencies are checked out on the following commits:

- `aocommon`: be6072490469977f513722b87f3032812438e968
- `eigen`: `3.3.7` (https://gitlab.com/libeigen/eigen/-/tags/3.3.7)
- `pybind11`: `v2.5.0` (https://github.com/pybind/pybind11/releases/tag/v2.5.0)

