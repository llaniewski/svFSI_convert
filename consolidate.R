library(reticulate)
library(dplyr)
library(Rcpp)

#  Import pyvista python module
pyvista = import('pyvista')

# LUMEN case
setwd("~/jijo/unstructured_mesh/VTK_lumen/")
surface_meshes_fn = c("LUMEN_INLET/lumen_inlet.vtp", "LUMEN_OUTLET/lumen_outlet.vtp", "LUMEN_WALL/lumen_wall.vtp")
volume_meshes_fn = c("lumen.vtu")

# WALL case
setwd("~/jijo/unstructured_mesh/VTK_wall/")
surface_meshes_fn = c("WALL_INLET/wall_inlet.vtp", "WALL_INNER/wall_inner.vtp", "WALL_OUTER/wall_outer.vtp", "WALL_OUTLET/wall_outlet.vtp")
volume_meshes_fn = c("wall.vtu")

region_id = 1
filenames = c(volume_meshes_fn, surface_meshes_fn)

# List of output files
outnames = c(
  "mesh-complete/mesh-complete.mesh.vtu",
  paste0("mesh-complete/mesh-surfaces/",basename(surface_meshes_fn))
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
cells = lapply( volume_sel, function(i) list(meshes[[i]]$n_cells, meshes[[i]]$cells, points$ref[(1+g_points_offset[i]):g_points_offset[i+1]]))
faces = lapply(surface_sel, function(i) list(meshes[[i]]$n_faces, meshes[[i]]$faces, points$ref[(1+g_points_offset[i]):g_points_offset[i+1]]))

# This implements an efficient method to find which faces (in the surfaces) are parts of which cells (of the volume mesh)
#  Rcpp is used, as otherwise this would be very slow
cppFunction('List find_cells(List faces, List cells) {
  List ret(cells.length() + faces.length());
  std::map<int, std::set<int> > points_in_cell;
  Rprintf("Constructing maps ...\\n");
  int celloffset=0;
  for (int m=0; m < cells.length(); m++) {
    List mesh = cells(m);
    int n_cells = mesh(0);
    IntegerVector tab = mesh(1);
    IntegerVector ref = mesh(2);
    IntegerVector cellids(n_cells);
    int cellidx = 0;
    for (int i=0; i < tab.length(); i += tab(i) + 1) {
      for (int j=0; j < tab(i); j++) {
        points_in_cell[ref(tab(i+j+1))].emplace(cellidx);
      }
      assert(cellidx < n_cells);
      cellids(cellidx) = cellidx + celloffset;
      cellidx++;
    }
    assert(cellidx == n_cells);
    ret(m) = cellids;
    celloffset += n_cells;
  }
  Rprintf("Inspecting maps ...\\n");
  for (int m=0; m < faces.length(); m++) {
    List mesh = faces(m);
    int n_faces = mesh(0);
    IntegerVector tab = mesh(1);
    IntegerVector ref = mesh(2);
    List cellids(n_faces);
    int faceidx = 0;
    for (int i=0; i < tab.length(); i += tab(i) + 1) {
      std::map<int, int> fc;
      for (int j=0; j < tab(i); j++) {
        std::set<int> c = points_in_cell[ref(tab(i+j+1))];
        for (auto k : c) fc[k]++;
      }
      std::vector<int> fi;
      for (auto k : fc) if (k.second == tab(i)) fi.push_back(k.first);
      assert(faceidx < n_faces);
      cellids(faceidx) = fi;
      faceidx++;
    }
    assert(cellidx == n_cells);
    ret(m + cells.length()) = cellids;
  }
  Rprintf("List ready\\n");
  return ret;
}')
ret = find_cells(faces, cells)

matches_one = sapply(surface_sel, function(i) all(sapply(ret[[i]],length) == 1))
if (!all(matches_one)) stop("Some surface faces matched more than one volume cell")

fc = do.call(c, lapply(ret, simplify2array)) + 1

# Assign the ids to the GlobalElementID and ModelRegionID fields in meshes
for (i in all_sel) {
  sel = g_cells_offset[i] + seq_len(n_cells[i])
  meshes[[i]]$cell_data$set_array(name="GlobalElementID", fc[sel])
  meshes[[i]]$cell_data$set_array(name="ModelRegionID", region_id)
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

