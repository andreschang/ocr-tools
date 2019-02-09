# OCR Tools
Open Climate Research is an ongoing project that aims to facilitate creative experimentation with modeled climate data. OCR Tools aims to be much more than a climate data viewer by enabling non-scientists to utilize a wide range of datasets and providing users with simple feedback conducive to learning. In addition to providing basic analysis functions, OCR Tools includes organizational and creative tools.

## Installing / Getting started

Run the following to install:
```python
pip install ocrtools
```

## Examples

- Open a NetCDF dataset with 

- ``` python
  import ocrtools as ocr
  cesm_TS = ocr.load('path/to/cesm_TS_data.nc', var='TS')
  ```

  If `var` is omitted, ocrtools will print out all variables in the dataset and ask you to specify a variable(s) of interest via command line. The dataset is then opened as an Xarray Dataset

- Create a `scope` object

  ``` python
  lima_peru = ocr.scope(location='Lima, Peru', yr0=1950, yrf=2000)
  ```

  * Location can also be specified by keyword arguments `lat_min`, `lat_max`, `lon_min`, and `lon_max`; or if none of these are given, location can be specified interactively by selecting areas on a map

- Subset your data

  ```python
  lima_TS = ocr.subset(cesm_TS, lima_peru)
  ```

- Select an area on a map and take the spatial average

  ```python
  from ocrtools import plt
  map_selection = ocr.scope()
  ```

  ```shell
  [OCR] Creating new scope object
  Enter yr0: 
  Enter yrf: 
  Select area(s) on map and close the pop-up window
  ```
  <img src="http://andreschang.com/unlinked/tk_selector_screenshot.png" width=70%/>

```shell
[OCR] Finished writing new scope object
```

```python
peru_TS = ocr.subset(cesm_TS, map_selection)
peru_avg_TS = ocr.spatial_average(peru_TS)
peru_avg_TS['TS'].plot()
plt.show()
```

