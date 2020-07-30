External dependencies are included via git submodules, see .gitmodules in the root of the source.

Please note that the dependencies are fixed to a specific commit, to make sure tests do not break as a 
result of silently updating the submodules. Dependencies are checked out on the following commits:

- `aocommon`: 05917d37e32a1bcd3fc985bc775b0b7b13d12532 (https://gitlab.com/aroffringa/aocommon)
- `eigen`: `3.3.7` (https://gitlab.com/libeigen/eigen/-/tags/3.3.7)
- `pybind11`: `v2.5.0` (https://github.com/pybind/pybind11/releases/tag/v2.5.0)

