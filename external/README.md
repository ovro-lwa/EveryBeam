External dependencies are included via git submodules, see .gitmodules in the root of the source.

Please note that the dependencies are fixed to a specific commit, to make sure tests do not break due to
a silent updateing of the submodules. Dependencies are checked out on the following commits:

- `aocommon`: 1f77c0c5e0d70507f227892353df8c5410cd87a4
- `eigen`: `3.3.7` (https://gitlab.com/libeigen/eigen/-/tags/3.3.7)
- `pybind11`: `v2.5.0` (https://github.com/pybind/pybind11/releases/tag/v2.5.0)

