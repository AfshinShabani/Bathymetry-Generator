Program Overview:
This Python-based program is designed to generate bathymetry files for the CE-QUAL-W2 hydrodynamic and water quality model, a widely used tool for simulating two-dimensional, longitudinal, and vertical water bodies like rivers, lakes, and estuaries. The program facilitates the creation of bathymetric data by processing topographic contour data and the centerline of a river or water body. By automating the bathymetry generation process, this tool helps modelers streamline the preparation of input data for CE-QUAL-W2.

Input Requirements:
To use this tool, users must provide three key input datasets:

1. Topographic Contours: The contours should be supplied as a shapefile, with an elevation attribute that the program will use to generate the bathymetry. It's important that the contour shapefile is free from geometry errors, such as dangling or self-intersecting lines, which could lead to processing issues.

2. River Centerline: This shapefile represents the main axis of the river. It must contain a 'Segment' column that correlates with the CE-QUAL-W2 segments and should not intersect with multiple W2 segments.

3. W2 Segments: A shapefile defining the W2 model segments, which must also include a 'Segment' column. Each segment should correspond to a CE-QUAL-W2 segment number, following the numbering conventions of the model (i.e., the first segment is 2, the second is 3, and so on).

Key Features and Assumptions:
The tool expects the W2 segments and centerline shapefiles to use CE-QUAL-W2's active segment numbering system. Please see the example file. 
Inactive cells are automatically accounted for by the program, meaning users do not need to include them in the input data—the program will insert inactive cells into the output where necessary.
This process ensures that the bathymetric data aligns correctly with the CE-QUAL-W2 model grid, providing an accurate representation of depth and riverbed structure for subsequent simulations.
By using this bathymetry generator, CE-QUAL-W2 modelers can save significant time in creating accurate bathymetric inputs and avoid potential errors in manually configuring bathymetric layers.

You can download the executable program for Windows and an example file here https://w2-bathymetry-generator.s3.us-east-2.amazonaws.com/W2+Bathymetry+Generator.zip
