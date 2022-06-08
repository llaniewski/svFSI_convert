# svFSI_convert
This script takes vtu volume mesh, and vtp surface meshes, matches points,
and cell ids and fills in GlobalNodeID, GlobalElementID, and ModelRegionID,
for it to work with svFSI

## Dependencies
In R:
```R
install.packages("reticulate")
install.packages("dplyr")
install.packages("Rcpp")
```

For python:
```shell
pip3 install pyvista
```
