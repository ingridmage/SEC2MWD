# SEC2MWD

A MATLAB toolbox for derivation of molecular weight distributions from size exclusion chromatography

# Description
Size exclusion chromatography (SEC) is a type of liquid chromatography used for separating molecules based on their size. This toolbox contains the pipeline for converting raw chromatograms to molecular weight distributions. The pipeline consists of seven steps, illustrated in the figure below.

The toolbox contains seven main functions, corresponding to each step in the calculation pipeline. Data is first imported into a MATLAB table where each row represents a SEC run and columns contain the raw data and additional meta data. For each subsequent step in the pipeline, the table is populated with more columns containing processed data and calculated results.  

In addition to the seven main functions, there is one supporting function that converts retention time values to molecular mass values. 

The modular architecture makes it easy to modify or replace individual steps of the pipeline.  


![image](https://github.com/ingridmage/SEC2MWD/assets/128492649/4194cc7e-276b-4e3f-a784-0bbbee6d6974)


# System requirements
SEC2MWD has been developed on a MATLAB R2022b, and earlier versions of MATLAB may not be supported. Only the base MATLAB installation is required, with no dependency on any of MATLAB's toolboxes.

# Use

The _samples_ folder contains an example script and data for running the workflow.

# License

The code is lisenced under GNU General Public License v3.0
