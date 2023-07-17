# TRSF (Topological Radio Source Finder)

[![Run tests](https://github.com/RhysAlfShaw/TRSF/actions/workflows/pytest.yaml/badge.svg)](https://github.com/RhysAlfShaw/TRSF/actions/workflows/pytest.yaml)
[![codecov](https://codecov.io/gh/RhysAlfShaw/TRSF/branch/main/graph/badge.svg?token=CFARAVMGZU)](https://codecov.io/gh/RhysAlfShaw/TRSF)

Currently in development tool for radio astronomers to find sources, specifically identifying components of well resolved sources.

!! CURRENTLY IN DEVELOPMENT !!

Source_props requires some compiled cython to use the faster region expansion method. To do this move into the source_props dir and use the following command.
```bash
python setup.py build_ext --inplace
```

Running in commandline:

```bash
python -m TRSF [Input Path] [cutout size] [sigma] [outputdir] [output file format]
```
