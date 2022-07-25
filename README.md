# SurfR
SurfR developing fork

## Gene2SProtein function

- **Description**: Detect Surface Proteins from a list of genes using surfy database
- **Required Packages**: import, rio, openxlsx
- **Inputs**:   
    1) vector of genes, 
    2) type of gene ID (e.g. hgnc, ensemble, entrez),
    3) if you want to print a tsv file 
    4) Surfy data-frame versioning (default newest, if not log) check size. 
- **Outputs**:  
    1) data.frame /(data.table) with filtered list of surface protein coding genes 
    2) .tsv file 
    3) standard output with number of surface protein found 
    4) Surfy data-frame used
- **Log**: 1) log with used Surfy data-frame 
- **Warnings**: 
    1) new Surfy data-frame version exist;  
    2) no surface protein found in your list of genes
- **Errors**: 1) no marched protein found in surfy database

