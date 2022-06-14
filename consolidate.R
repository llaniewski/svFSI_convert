library(reticulate)
library(dplyr)
library(Rcpp)

#  Set working directory:
#setwd("~/jijo/svFSI_mesh_conversion_jijo/mesh_files_converted_with_openFOAM/")

#  Import pyvista python module
pyvista = import('pyvista')

#  List of input files
surface_meshes_fn = c("LUMEN_INLET/LUMEN_INLET.vtp", "LUMEN_OUTLET/LUMEN_OUTLET.vtp", "LUMEN_WALL/LUMEN_WALL.vtp")
volume_meshes_fn = c("lumen.vtu")


filenames = c(volume_meshes_fn, surface_meshes_fn)
outnames = paste0("out/",filenames)

# List of output files
outnames = c(
  "mesh-complete/mesh-complete.mesh.vtu",
  "mesh-complete/mesh-surfaces/lumen_inlet.vtp",
  "mesh-complete/mesh-surfaces/lumen_outlet.vtp",
  "mesh-complete/mesh-surfaces/lumen_wall.vtp"
)

# Read meshes
meshes = lapply(filenames, pyvista$read)

# Indexes of the volume and surface meshes
volume_sel  = seq_along( volume_meshes_fn)
surface_sel = seq_along(surface_meshes_fn) + length(volume_meshes_fn)
all_sel = seq_along(meshes)

# Some helper vectors
n_points = sapply(meshes, function(m) m$n_points)
g_points_offset = c(0,cumsum(n_points))
n_cells = sapply(meshes, function(m) m$n_cells)
g_cells_offset = c(0,cumsum(n_cells))

# This part constructs a big table of points and assignes them unique ids
points = lapply(meshes, function(m) m$points)
points = do.call(rbind,points)
points = data.frame(points)
coords = names(points)
points$id = 1:nrow(points)
points = points %>% group_by_at(coords) %>% mutate(ref = first(id))
refs = unique(points$ref)
nid = rep(NA, nrow(points))
nid[refs] = seq_along(refs)
points$ref = nid[points$ref]

# Assignes these ids (in "ref" column) to GlobalNodeID field in meshes
for (i in all_sel) {
  meshes[[i]]$point_data$set_array(name="GlobalNodeID", points$ref[g_points_offset[i] + seq_len(n_points[i])])
}

# Now the hardest part

do.rbind.fillzero = function(x) {
  n = max(sapply(x,ncol))
  do.call(rbind, lapply(x, function(tab) {
    if (ncol(tab) < n) {
      cbind(tab, matrix(0,nrow(tab),n-ncol(tab)))
    } else {
      tab
    }
  }))
}

# This part gets cells from volume meshes and faces from surface meshes as arrays
#  the first column of the arrays are the number of points in respective cells/faces
#  eg. [3 10 11 12 0 0] means a 3 point cell (triangle) with vertex indexes 10, 11, 12
cells = lapply(volume_sel,function(i) {
  ind = meshes[[i]]$cells
  pad = length(ind) / meshes[[i]]$n_cells
  ind = matrix(ind,ncol=pad,byrow = TRUE)
  for (j in 1:max(ind[,1])) {
    sel = j <= ind[,1]
    ind[sel,j+1] = points$ref[ind[sel,j+1]+1+g_points_offset[i]]
  }
  ind
})
cells = do.rbind.fillzero(cells)

faces = lapply(surface_sel,function(i) {
  ind = meshes[[i]]$faces
  pad = length(ind)/meshes[[i]]$n_faces
  ind = matrix(ind,ncol=pad,byrow = TRUE)
  for (j in 1:max(ind[,1])) {
    sel = j <= ind[,1]
    ind[sel,j+1] = points$ref[ind[sel,j+1]+1+g_points_offset[i]]
  }
  ind
})
faces = do.rbind.fillzero(faces)


# This implements an efficient method to find which faces (in the surfaces) are parts of which cells (of the volume mesh)
#  Rcpp is used, as otherwise this would be very slow
cppFunction('List find_cells(IntegerMatrix faces, IntegerMatrix cells) {
  List ret(faces.nrow());
  std::map<int, std::set<int> > points_in_cell;
  Rprintf("Constructing maps ...\\n");
  for (int i=0; i<cells.nrow(); i++) {
    for (int j=0; j<cells(i,0); j++) {
      points_in_cell[cells(i,1+j)].emplace(i);
    }
  }
  Rprintf("Inspecting maps ...\\n");
  for (int i=0; i<faces.nrow(); i++) {
    std::map<int, int> fc;
    for (int j=0; j<faces(i,0); j++) {
      std::set<int> c = points_in_cell[faces(i,1+j)];
      for (auto k : c) fc[k]++;
    }
    std::vector<int> fi;
    for (auto k : fc) if (k.second == faces(i,0)) fi.push_back(k.first);
    ret[i] = fi;
  }
  Rprintf("List ready\\n");
  return ret;
}')
ret = find_cells(faces, cells)

if (any(sapply(ret,length)!= 1)) stop("Some surface faces matched more than one volume cell")

fc = simplify2array(ret) + 1  # +1 to go from C indexing to R/fortran indexing

# Make a list with ids and region id's for the cells and faces
elements = data.frame(id = 1:nrow(cells), region = rep(volume_sel, times=n_cells[volume_sel]))
elements = rbind(elements,elements[fc,])
if (nrow(elements) != nrow(cells)+nrow(faces)) stop("Something went wrong")

# Assign the ids to the GlobalElementID and ModelRegionID fields in meshes
for (i in all_sel) {
  sel = g_cells_offset[i] + seq_len(n_cells[i])
  meshes[[i]]$cell_data$set_array(name="GlobalElementID", elements$id[sel])
  meshes[[i]]$cell_data$set_array(name="ModelRegionID", elements$region[sel])
}

# Writing the meshes back to the files
for (i in all_sel) {
  fn = outnames[i]
  p = dirname(fn)
  if (!dir.exists(p)) dir.create(p,recursive = TRUE)
  cat("Writing",fn,"\n")
  # meshes[[i]]$save(fn)  # this would be normally do, but does not support "append" mode
  ext = paste0(".", tools::file_ext(fn))
  writer = meshes[[i]]$'_WRITERS'[[ext]]()
  writer$SetDataModeToAppended()
  writer$SetFileName(fn)
  writer$SetInputData(meshes[[i]])
  writer$Write()
  rm(writer)
}

