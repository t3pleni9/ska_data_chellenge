# Setup
* Install conda
* For Linux based machine

  `conda env create -f environment.yml`
* For Mac based machines
  * install [sofia-2](https://github.com/SoFiA-Admin/SoFiA-2)
  
  `conda env create -f mac-environment.yml`
  
  
# Merging Catalogues

`$ PYTHONPATH=. ipython`


```from skavengers.catalogue import Catalogue 
    c = Catalogue('~/Downloads/sky_dev_truthcat_v1.1.txt',  'ascii')
    c1 = Catalogue(path='./output/sky_fits_sources_1_cat.sql', input_type='sql')
    c2 = Catalogue(path='./output/sky_fits_sources_2_cat.sql', input_type='sql')
    c3 = Catalogue(path='./output/sky_fits_sources_3_cat.sql', input_type='sql')
    c1.merge_catalogues([c2, c3], 1)```


  
  
