from posydon.grids.psygrid import PSyGrid

grid_dir = "/home/IITB/CompnalAstrphyNRelvity/23n0315/MESA/mesa-24.08.1/binary_pop/runs"
output_file = "/home/IITB/CompnalAstrphyNRelvity/23n0315/posydon/POSYDON_2/POSYDON/posydon/grids/binary_grid/my_grids"

mygrid = PSyGrid(verbose=True).create(grid_dir, output_file)
