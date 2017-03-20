"""Obtain mutation data from TCGA

Raw data files were obtained from https://gdc-portal.nci.nih.gov/search/f?filters=%7B%22op%22:%22and%22,%22content%22:%5B%7B%22op%22:%22%3C%3D%22,%22content%22:%7B%22field%22:%22cases.diagnoses.age_at_diagnosis%22,%22value%22:%5B7305%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.data_category%22,%22value%22:%5B%22Simple%20Nucleotide%20Variation%22%5D%7D%7D,%7B%22op%22:%22in%22,%22content%22:%7B%22field%22:%22files.access%22,%22value%22:%5B%22open%22%5D%7D%7D%5D%7D&pagination=%7B%22files%22:%7B%22from%22:0,%22size%22:20,%22sort%22:%22data_category:desc,%22%7D%7D&facetTab=files
Data include all publicly-available simple nucleotide variation information on the website. Descriptions of raw data
files are available at https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification


"""