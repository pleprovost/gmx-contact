# Contact 
Calucate the pair distance between two selection and count as contact distances within cutoff.

## Requirement
* Gromacs 2019 or more
* C++17

## Compile
Copy the cmake directory from /whatever/gromacs/share/gromacs/cmake 
```
mkdir build
cd build
cmake ..
make 
```
## Usage
```
/path/to/contact -f traj.xtc -s prod.tpr -o -tu ns
```

