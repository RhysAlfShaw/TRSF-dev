# TRSF (Topological Radio Source Finder)

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
