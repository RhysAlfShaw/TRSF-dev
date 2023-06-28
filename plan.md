# Plan for TRSF Python Package

This is a development plan for creating TRSF project.

1. 2D continum Image (currently GUDHI implementation)

2. Source subtraction functionality.

3. Multigridding functionality.

4. Cubical Ripser Persistent Homology implementation (for speed up and more convenient data structure).

5. Initial Data Cube Source finding Implementation.

6. Option to generate a source catalogue (saving to fits, hdf or csv file).

## Input Requirements.

- file_name and location
- Area of interest in image (in pixels or sky postion).
- return a component Catalogue.
- options for data return.

## Output Style.

- Finds point sources and others as gaussians (ellipse gaussian parameters for fitting).
- Extended emision that encompases other point sources will be flagged as such (will this be fitted with a ellipse around the polygon or other).
- Intergrated flux.
- Residual from subtraction, cutout or whole image.

## Processing Pipeline

- create subimages if necessary.
- find source with cubical complex's.
- fit all small sources with a 2d gaussian, after regions are maximized to allow for better fit outcomes.
- larger diffuse regions remain as polygons surrounding them.

## In future a system for make associations or flagging spurious sources.

## Enhancements.

- Error on gaussian fitting.
- Optimazation with C code where appropriate and if necessary in the future.
