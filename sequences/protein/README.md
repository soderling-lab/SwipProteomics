# SWIP's predicted structure

A structure of the SWIP protein was generated using [Phyre](http://www.sbg.bio.ic.ac.uk/phyre2/html/page.cgi?id=index).
The SWIP protein has no recognizable domains, but is predcited to be highly helical. 
It's structure appears to be modular, organized into two major domains. 
SWIP interacts directly with Strumpellin, this interaction is thought to be
mediated through its C-terminus ([Jia et al., 2010](../refs/Jia_2010)). The
SWIP<sup>P1019<\sup> mutation is found in the C-terminal domain and is predicted
to be structurally disruptive [(Missense3d)](http://www.sbg.bio.ic.ac.uk/~missense3d/). 

## Creating a gif of the SWIP protein
In pymol, generate a gif of the SWIP protein can be generated using the following commands:
```
# Generate a gif.
mset 1,180 
util.mroll 1, 180, 1
set ray_trace_frames, 1 
set cache_frames, 0

# Save as a video: File > save as video > PNG images
# The image stack was made into a gif with the img2gif shell script.

```

The full use of Pymol without a liscense is restricted.
A better alternative to pymol is chimera. Convert `protein.pse` (pymol specific format) 
to the univeral .pdb format and create a video in chimera.

```
# Convert the video to a stack  of images.
ffmpeg -i infile.avi -f image2 image-%03d.jpg

# Convert the images to a gif.
img2gif JPG/ 1 -1 wt.gif
```
