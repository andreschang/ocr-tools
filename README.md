# OCR Tools
Open Climate Research is an ongoing project that aims to facilitate creative experimentation with modeled climate data. OCR Tools aims to be much more than a climate data viewer by enabling non-scientists to utilize a wide range of datasets and providing users with simple feedback conducive to learning. In addition to providing basic analysis functions, OCR Tools includes organizational and creative tools.

## Installing / Getting started

Make sure your Python environment meets a few basic requirements. OCR Tools just needs Numpy, Matplotlib, Basemap, and NetCDF4. 

[Anaconda](https://conda.io/docs/user-guide/tasks/manage-environments.html) users can set up their environment like this:

```shell
conda env create -f environment.yml
source activate ocrtools
```

Next, set your data directories in the Global Variables section of reformat_data.py (around line 20), and you're all set. Reference the reformat_data.py from any python file in the same directory to use tools. A few sample usages can be found in reformat_live.py.

## Features

* Explore a climate data netCDF without knowing anything about climate data or netCDFs.
* Perform analysis with just a few basic parameters – at minimum, `variable`, `yr0`, and `yrf`.
* Plot a timeseries of the spatial average of any variable over a user-defined area and print a map of the area.
* Reformat giant netCDFs into a set of manageable chunks that can be reassembled automatically by OCR Tools when needed.
* Automatic file naming and organization

## Contributing

Contributors are more than welcome! There is a lot of work to be done, and ultimately, I'd like OCR Tools to be fully functional from within a GUI, so that users don't need to interface with any code at all (if they don't want to). Please reach out to me at andresdanielchang@gmail.com if you're interested.