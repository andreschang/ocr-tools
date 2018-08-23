# OCR Tools
Open Climate Research is an ongoing project that aims to facilitate creative experimentation with modeled climate data. OCR Tools aims to be much more than a climate data viewer by enabling non-scientists to utilize a wide range of datasets and providing users with simple feedback conducive to learning. In addition to providing basic analysis functions, OCR Tools includes organizational and creative tools.

## Installing / Getting started

OCR Tools runs in Python, and it requires a few basic packages: [Matplotlib](https://matplotlib.org/), [Basemap](https://matplotlib.org/basemap/), and [NetCDF4](http://unidata.github.io/netcdf4-python/). The packages can be installed individually using the [Anaconda Navigator](https://www.anaconda.com/distribution/) or [Canopy](https://www.enthought.com/product/canopy/) desktop applications (no command line needed) . 

Conda users may also set up their environment through the command line by typing:

```shell
conda env create -f environment.yml
```

and then switching to the new environment by typing ```source activate ocrtools``` (Mac) or ```activate ocrtools``` (Windows).

Once your environment is ready, start exploring ```demo.py```.  The ```ocrtools``` folder is preloaded with a ```data``` directory and sample climate data, so it should run out of the box.



## Features

* Explore a climate data netCDF without knowing anything about climate data or netCDFs.
* Perform analysis with just a few basic parameters – at minimum, `variable`, `yr0`, ```src```, and `yrf`.
* Plot a timeseries of the spatial average of any variable over a user-defined area and print a map of the area.
* Reformat giant netCDFs into a set of manageable chunks that can be reassembled automatically by OCR Tools when needed.
* Automatic file naming and organization

## Contributing

Contributors are more than welcome! There is a lot of work to be done, and ultimately, I'd like OCR Tools to be fully functional from within a GUI, so that users don't need to interface with any code at all (if they don't want to). Please reach out to me at andresdanielchang@gmail.com if you're interested.