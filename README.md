# GJ3929b-MIRI-F1500W-Reproducibility
Steps to reproduce the results from Connors et al 2026 using Eureka! and Erebus for the Rocky Worlds DDT Data Challenge:

Download the contents of this repo.

Download the data from [Kaggle](https://www.kaggle.com/datasets/stsci/rocky-worlds-gj-3929b-observations) and unzip it into the root directory of this repo.

Erebus is hard coded to treat each individual visit as its own eclipse light curve based on the unique visit ID. To get around this for the split observation of eclipse 4, rename the following files:
`jw09235005001_02101_00001-seg001_mirimage_uncal.fits` to `jw09235004001_03101_00001-seg005_mirimage_uncal.fits`
`jw09235005001_02101_00001-seg002_mirimage_uncal.fits` to `jw09235004001_03101_00001-seg006_mirimage_uncal.fits`
`jw09235005001_02101_00001-seg003_mirimage_uncal.fits` to `jw09235004001_03101_00001-seg007_mirimage_uncal.fits`

Erebus is also hardcoded to expect the eclipses to be in a folder named `JWST`. In the `GJ_3929b_Observations` create a subfolder named `JWST` and move the 4 eclipse folders into it (such that the directory structure is):

```
GJ3929b-MIRI-F1500w-Reproducibility
  > run_gj3929b.py
  > gj3929b_cfg.yaml
  > gj3929b_planet.yaml
  > GJ_3929b_Observations
    > JWST
      > eclipse1
      > eclipse2
      > eclipse3
      > eclipse4
```

Open up `run_gj3929b.py` and update the following lines at the top of the file
```python
crds_path = #'/home/nconnors/crds_cache'
uncal_path = #'/home/nconnors/Research/GJ3929b_analysis/GJ_3929b_Observations'
```
Put in the absolute path to your crds cache (required to run Eureak! to process Stage 1 and 2).
Put in the absolute path to your GJ_3929b_Observations folder.

Now create a conda environment with python 3.12 and install Erebus version 0.8.1.

```
conda create --name erebus_env
conda activate erebus_env
conda install python=3.12
pip install erebus-exoplanet==0.8.1
```

Also install Eureka! following the instructions [here](https://eurekadocs.readthedocs.io/en/stable/installation.html). The exact version used in the Kaggle submission is v1.2.2.

Now simply call `python run_gj3929b.py`.

# TODO: Show how to get posterior distributions for results.
