# Contact 
Calucate the pair distance between two selection and count as contact distance within cutoff

# Compile
Copy the cmake directory from /whatever/gromacs/share/gromacs/cmake 

mkdir build
cd build
cmake ..
make 

# Usage
/path/to/contact -f traj.xtc -s prod.tpr -o -tu ns

