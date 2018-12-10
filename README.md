# tsdf
### Overview
Code to fuse multiple Kinect point clouds into a truncated signed distance function (TSDF) voxel volume. Visualization of 3D surface meshes with Marching Cubes algorithm. 
### Usage: 
`
./tsdf [options]`

* `-p:                 point cloud prefix`
* `-m  <int>:          maximum number of scans that should be processed`
* `-e  <int>:          process only every e-th scan`
* `-r  <float>:        resolution of TSDF map in meters (default 0.01)`
* `-td <float>:        trunctation distance in meters (default 0.02)`
* `-minw <float>:      minimum weight for marching cubes (default 0.0)`
* `-o <output.ply>:    output mesh file name (default mesh.ply)`
* `-c:                 integrate color (default false)`
* `-csv :              write csv file (default false)`
* `-ocsv <output.csv>: output csv file name (default tsdf.csv)`
* `-h:                 help`
`
#### Example: 
```
./tsdf -p ../data/fr1-room/ncloud -m 10 -r 0.01 -td 0.02 -c

```
