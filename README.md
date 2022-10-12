Code for the manuscript  
## **Global determinants of insect genetic diversity** <br/>

Where to find things:

## Data
*Need to update once published. Some files are too large to host on github and will need to be kept in a remote server that will be provided upon acceptance*  


Raw climate data:  
`data/climate_raw`  

Raster templates for aggregating data:  
`data/templates`

Climate data aggregated for analysis:  
`data/climate_agg`

Non-raster climate data:  
`data/climate_poly`

Invasive species lists:  
`data/invasive`

## Metadata
Metadata with information on predictor variables and invasive species lists is provided in the `meta/` folder.  

## Notebooks
These are analysis notebooks containing annotated code.  `.ipynb` are jupyter notebooks (open using python jupyter notebook), `.Rmd` are RMarkdown documents, and `.qmd` are Quarto documents (open `.Rmd` and `.qmd` in RStudio).  

Aggregating enviromental data for analysis:  
`step-0_aggregate-rasters.ipynb`

Filter and align sequence data:  
`step-1_seq-filter-align-test.ipynb`  

Calculate statistics from sequence data:  
`step-2_seq-stats-test.ipynb`  

Exploring, wrangling, and filtering data before modeling:  
`step-3_data-exploration.Rmd`

Modeling:  
`step-4-models.qmd`  

Publication figures:  
`publication_figures.Rmd`

These are the reports necessary to regenerate the analyses for publication. All other `.Rmd` files contain code for data exploration that are unnecessary to recreate the analysis.  

## Scripts
All helper R scrips are contained in the `R/` folder. They are called in their respective report.  


## Output
All files in the `output/` folder and subfolders are able to be regenerated with the code provided.  
