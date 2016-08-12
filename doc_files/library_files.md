###sasmol files

The library contains files with classes and methods to define molecular objects, read and write PDB and binary dcd files, calculate molecular properties, create and merge subsets from base molecules.  

Files with methods are provided to accommodate CHARMM force-field atom naming, dictionaries with atomic properties including x-ray and neutron scattering and volumetric values.

Utility methods to handle vector and matrix operations, overlap transformations, and a method to interactively transfer coordinates from python to VMD for visualization are provided.


---

###[sasmol.py](sasmol/sasmol_overview.md) 
 
 - contains the SasMol class to define sasmol objects (setters and getters)

###[file_io.py](file_io/file_io_overview.md)

 - contains the Files class for PDB/DCD input / output

###[calculate.py](calculate/calculate_overview.md)

 - contains methods to calculate properties from coordinates

###[operate.py](operate/operate_overview.md)

 - contains methods to perform geometric transformations on sasmol objects

###[subset.py](subset/subset_overview.md)

 - contains methods to extract / combine sub-molecules from / to other sasmol instances

###[linear_algebra.py](linear_algebra/linear_algebra.md)

 - contains mathematical methods for handling molecular coordinates
 

###[topology.py](topology/topology_overview.py)

 - contains methods to accommodate force-field naming conventions

###[properties.py](properties/properties_overview.md)

 - contains dictionaries of atomic properties


###[view.py](view/view_overview.md)

 - contains convience methods to link coordinates to / from VMD
