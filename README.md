# Code for the manuscript: **Global determinants of insect genetic diversity**

Currently on bioRxiv:
> **Global determinants of insect genetic diversity**\
Connor M French, Laura D Bertola, Ana C Carnaval, Evan P Economo, Jamie M Kass, David J Lohman, Katharine A Marske, Rudolf Meier, Isaac Overcast, Andrew J. Rominger, Phillip Staniczenko, Michael J Hickerson *bioRxiv* 2022.02.09.479762; doi: <https://doi.org/10.1101/2022.02.09.479762>

Where to find things:

## Data

All climate data is available online at the links provided in Supplementary Table 1 or `meta/raster_list.csv`. The BOLD insect data is available at **figshare_link** as `bold_data_insects.txt`, which belongs in the `data/` folder.

The only data subfolders that come pre-populated are `data/templates`, which contains raster templates for aggregating raw climate data, and `data/invasive`, which contains spreadsheets with invasive species information for filtering out invasive species from the data set.

For the scripts to run, you need to create the following folders:

For raw climate `data/climate_raw`

Climate data aggregated for analysis:\
`data/climate_agg`

Non-raster climate `data/climate_poly`

Folder to put output:\
`output/`

-   all necessary output subfolders will be generated by code

## Metadata

Metadata with information on predictor variables and invasive species lists is provided in the `meta/` folder.

## Notebooks

These are analysis notebooks containing annotated code. `.ipynb` are jupyter notebooks (open using python jupyter notebook), `.Rmd` are RMarkdown documents, and `.qmd` are Quarto documents (open `.Rmd` and `.qmd` in RStudio or VS Code). The links go to the readable jupyter notebook file or github markdown file so you don't have to download the repo to view the code. Below each link is the path to the raw report file.

Aggregating environmental data for analysis:\
[step-0_aggregate-rasters](reports/step-0_aggregate-rasters.ipynb)

`reports/step-0_aggregate-rasters.ipynb`

Filter and align sequence:   [step-1_seq-filter-align-test](reports/step-1_seq-filter-align-test.ipynb)

`reports/step-1_seq-filter-align-test.ipynb`

Calculate statistics from sequence:  
[step-2_seq-stats-test](reports/step-2_seq-stats-test.ipynb)

`reports/step-2_seq-stats-test.ipynb`

Exploring, wrangling, and filtering data before modeling:\
[step-3_data-exploration](reports/step-3_data-exploration.md)

`reports/step-3_data-exploration.Rmd`

Modeling:\
[step-4-models](reports/step-4-models.md)

`reports/step-4-models.qmd`

Publication figures:\
[publication_figures](reports/publication_figures.md)

`reports/publication_figures.Rmd`

-   I didn't optimize the plotting output for the report, so sorry for some unreadable plots

These are the reports necessary to regenerate the analyses for publication.

## Scripts

All helper R scrips are contained in the `R/` folder. They are called in their respective report.
