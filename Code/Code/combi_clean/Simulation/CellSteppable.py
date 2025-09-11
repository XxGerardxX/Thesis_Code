from cc3d.core.PySteppables import *
from fipy import CellVariable, Grid3D
import all_values as al
import sys
from collections import defaultdict
global mesh
from numpy import *
import numpy as np
from builtins import range
import random

global cytokines
global cellpresent
global fullFileName
# globals for search
global setLambda
global saturation_coefficient




sys.modules['__main__'] = sys.modules[__name__]

# Predefinition of parameters


setlambda = 2000
saturation_coefficient = 10 ** -6


def randopos_3D(typein, layer, xboundary=None, yboundary=None, zboundary=None,
                lim=None, cell_freeze=None, nx=al.nx, ny=al.ny, nz=al.nz):
    """
    3D version of randomization of position within lattice for burn injury simulation
    
    Args:
        typein: Placement type (0=systemic, 1=wound edges, 2=local wound, 3=systemic-wound)
        xboundary: X boundary parameter (default: nx//4)
        yboundary: Y boundary parameter (default: ny//4) 
        zboundary: Z boundary parameter for tissue layers (default: nz//4)
        lim: Limit parameter for boundary zones (default: 10)
        cell_freeze: Cell freeze parameter (default: 5)
        nx, ny, nz: Field dimensions
        
    Returns:
        list: [x, y, z] coordinates
    """
    
    # Set default parameters if not provided
    if xboundary is None:
        xboundary = nx // 4
    if yboundary is None:
        yboundary = ny // 4
    if zboundary is None:
        zboundary = nz // 4  # Tissue layer boundary
    if lim is None:
        lim = 10
    if cell_freeze is None:
        cell_freeze = 2
    
    if typein == 0:
        # Systemic placement (boundary regions in 3D space)
        randxhi = random.uniform(0, nx)
        randxlo = random.uniform(0, nx)
        randyall = random.uniform(0, ny)
        randyhi = random.uniform(0, ny)
        randylo = random.uniform(0, ny)
        randxall = random.uniform(0, nx)
        
        # Add z-dimension variations for systemic placement
        randzhi = random.uniform(0, nz)
        randzlo = random.uniform(0, nz)
        randzall = random.uniform(0, nz)
        randzmid = random.uniform(zboundary, nz - zboundary)
        
        # Create 3D placement options (edges and faces of the 3D space)
        hixvals = [randxhi, randyall, randzall]
        loxvals = [randxlo, randyall, randzall]
        hiyvals = [randxall, randyhi, randzall]
        loyvals = [randxall, randylo, randzall]
        hizvals = [randxall, randyall, randzhi]
        lozvals = [randxall, randyall, randzlo]
        midvals = [randxall, randyall, randzmid]  # Mid-tissue placement
        
        rands = [hixvals, loxvals, hiyvals, loyvals, hizvals, lozvals, midvals]
        
        randselect = random.randint(0, len(rands) - 1)
        randx = int(rands[randselect][0])
        randy = int(rands[randselect][1])
        randz = int(rands[randselect][2])
        
        # Clamp values to prevent out of bounds
        randx = max(0, min(randx, nx - 1))
        randy = max(0, min(randy, ny - 1))
        randz = max(0, min(randz, nz - 1))
    
    elif typein == 1:
        # Wound edges placement (3D boundary regions around wound)
        # Left boundary region
        randxhi1 = random.uniform(max(0, xboundary - lim - cell_freeze), min(nx, xboundary + lim + cell_freeze))
        randxlo1 = random.uniform(max(0, xboundary - lim), min(nx, xboundary + lim + cell_freeze))
        randyall1 = random.uniform(max(0, yboundary + lim), min(ny, ny - lim - yboundary))
        
        # Lower boundary region
        randyhi1 = random.uniform(max(0, yboundary - lim - cell_freeze), min(ny, yboundary + lim + cell_freeze))
        randylo1 = random.uniform(max(0, yboundary - lim - cell_freeze), min(ny, yboundary + lim + cell_freeze))
        randxall1 = random.uniform(max(0, xboundary + lim), min(nx, nx - lim - xboundary))
        
        # Right boundary region
        randxhi2 = random.uniform(max(0, nx - xboundary - cell_freeze), min(nx, nx - xboundary - lim + cell_freeze))
        randxlo2 = random.uniform(max(0, nx - xboundary - cell_freeze - lim), min(nx, nx - xboundary))
        randyall2 = random.uniform(max(0, yboundary + lim), min(ny, ny - lim - yboundary))
        
        # Upper boundary region
        randyhi2 = random.uniform(max(0, ny - yboundary - lim - cell_freeze), min(ny, ny - yboundary + lim + cell_freeze))
        randylo2 = random.uniform(max(0, ny - yboundary - cell_freeze - lim), min(ny, ny - yboundary + lim))
        randxall2 = random.uniform(max(0, xboundary + lim), min(nx, nx - lim - xboundary))
        
        # Z-dimension wound edge variations (tissue layer boundaries)
        # Superficial layer edges (epidermis-dermis boundary)
        randzsurf = random.uniform(max(0, zboundary - lim - cell_freeze), min(nz, zboundary + lim + cell_freeze))
        # Deep layer edges (dermis-subcutaneous boundary)  
        randzdeep = random.uniform(max(0, nz - zboundary - lim - cell_freeze), min(nz, nz - zboundary + lim + cell_freeze))
        # Mid-layer placement
        randzmid = random.uniform(max(0, zboundary + lim), min(nz, nz - zboundary - lim))



        # 3D wound edge positions
        hixvals1 = [randxhi1, randyall1, randzsurf]
        loxvals1 = [randxlo1, randyall1, randzsurf]
        hiyvals1 = [randxall1, randyhi1, randzsurf]
        loyvals1 = [randxall1, randylo1, randzsurf]
        hixvals2 = [randxhi2, randyall2, randzdeep]
        loxvals2 = [randxlo2, randyall2, randzdeep]
        hiyvals2 = [randxall2, randyhi2, randzdeep]
        loyvals2 = [randxall2, randylo2, randzdeep]
        
        # Z-boundary specific placements
        zsurf_vals = [random.uniform(max(0, xboundary), min(nx, nx - xboundary)), 
                     random.uniform(max(0, yboundary), min(ny, ny - yboundary)), randzsurf]
        zdeep_vals = [random.uniform(max(0, xboundary), min(nx, nx - xboundary)), 
                     random.uniform(max(0, yboundary), min(ny, ny - yboundary)), randzdeep]
        zmid_vals = [random.uniform(max(0, xboundary), min(nx, nx - xboundary)), 
                    random.uniform(max(0, yboundary), min(ny, ny - yboundary)), randzmid]
        
        rands = [hixvals1, loxvals1, hiyvals1, loyvals1, hixvals2, loxvals2, hiyvals2, loyvals2,
                zsurf_vals, zdeep_vals, zmid_vals]
        
        randselect = random.randint(0, len(rands) - 1)
        randx = int(rands[randselect][0])
        randy = int(rands[randselect][1])
        randz = int(rands[randselect][2])




    
    elif typein == 2:
        if layer == 0:
            # Local wound placement (interior wound region in 3D)
            randxhi = random.uniform(max(0, xboundary + lim + cell_freeze), min(nx, nx - lim - xboundary - cell_freeze))
            randxlo = random.uniform(max(0, xboundary + lim + cell_freeze), min(nx, nx - lim - xboundary - cell_freeze))

            randyall = random.uniform(max(0, yboundary + lim), min(ny, ny - lim - yboundary))
            randyhi = random.uniform(max(0, yboundary + lim + cell_freeze), min(ny, ny - lim - yboundary - cell_freeze))
            randylo = random.uniform(max(0, yboundary + lim + cell_freeze), min(ny, ny - lim - yboundary - cell_freeze))
            randxall = random.uniform(max(0, xboundary + lim), min(nx, nx - lim - xboundary))

            # Z-dimension for local wound (focus on specific tissue layers)
            # Epidermis layer (0-20% of tissue depth)
            randzepi = random.uniform(0, max(1, zboundary))
            # Dermis layer (20-70% of tissue depth)
            randzder = random.uniform(max(0, zboundary), min(nz, nz - zboundary))
            # Subcutaneous layer (70-100% of tissue depth)
            randzsub = random.uniform(max(0, nz - zboundary), max(1, nz - cell_freeze))

            # Local wound 3D positions across different tissue layers
            hixvals_epi = [randxhi, randyall, randzepi]
            loxvals_epi = [randxlo, randyall, randzepi]
            hiyvals_epi = [randxall, randyhi, randzepi]
            loyvals_epi = [randxall, randylo, randzepi]

            hixvals_der = [randxhi, randyall, randzder]
            loxvals_der = [randxlo, randyall, randzder]
            hiyvals_der = [randxall, randyhi, randzder]
            loyvals_der = [randxall, randylo, randzder]

            hixvals_sub = [randxhi, randyall, randzsub]
            loxvals_sub = [randxlo, randyall, randzsub]
            hiyvals_sub = [randxall, randyhi, randzsub]
            loyvals_sub = [randxall, randylo, randzsub]

            rands = [hixvals_epi, loxvals_epi, hiyvals_epi, loyvals_epi,
                    hixvals_der, loxvals_der, hiyvals_der, loyvals_der,
                    hixvals_sub, loxvals_sub, hiyvals_sub, loyvals_sub]

            randselect = random.randint(0, len(rands) - 1)
            randx = int(rands[randselect][0])
            randy = int(rands[randselect][1])
            randz = int(rands[randselect][2])

        if layer == 1:
            # Local wound placement (interior wound region in 3D)
            randxhi = random.uniform(max(0, xboundary + lim + cell_freeze), min(nx, nx - lim - xboundary - cell_freeze))
            randxlo = random.uniform(max(0, xboundary + lim + cell_freeze), min(nx, nx - lim - xboundary - cell_freeze))
            randyall = random.uniform(max(0, yboundary + lim), min(ny, ny - lim - yboundary))
            randyhi = random.uniform(max(0, yboundary + lim + cell_freeze), min(ny, ny - lim - yboundary - cell_freeze))
            randylo = random.uniform(max(0, yboundary + lim + cell_freeze), min(ny, ny - lim - yboundary - cell_freeze))
            randxall = random.uniform(max(0, xboundary + lim), min(nx, nx - lim - xboundary))

            # Z-dimension for local wound (focus on specific tissue layers)
            # Epidermis layer (0-20% of tissue depth)
            randzepi = random.uniform(0, max(1, zboundary))
            # Local wound 3D positions across different tissue layers
            hixvals_epi = [randxhi, randyall, randzepi]
            loxvals_epi = [randxlo, randyall, randzepi]
            hiyvals_epi = [randxall, randyhi, randzepi]
            loyvals_epi = [randxall, randylo, randzepi]

            rands = [hixvals_epi, loxvals_epi, hiyvals_epi, loyvals_epi]

            randselect = random.randint(0, len(rands) - 1)
            randx = int(rands[randselect][0])
            randy = int(rands[randselect][1])
            randz = int(rands[randselect][2])

        elif layer == 2:
            # Local wound placement (interior wound region in 3D)
            randxhi = random.uniform(max(0, xboundary + lim + cell_freeze), min(nx, nx - lim - xboundary - cell_freeze))
            randxlo = random.uniform(max(0, xboundary + lim + cell_freeze), min(nx, nx - lim - xboundary - cell_freeze))
            randyall = random.uniform(max(0, yboundary + lim), min(ny, ny - lim - yboundary))
            randyhi = random.uniform(max(0, yboundary + lim + cell_freeze), min(ny, ny - lim - yboundary - cell_freeze))
            randylo = random.uniform(max(0, yboundary + lim + cell_freeze), min(ny, ny - lim - yboundary - cell_freeze))
            randxall = random.uniform(max(0, xboundary + lim), min(nx, nx - lim - xboundary))

            # Z-dimension for local wound (focus on specific tissue layers)
            # Dermis layer (20-70% of tissue depth)
            randzder = random.uniform(max(0, zboundary), min(nz, nz - zboundary))

            hixvals_der = [randxhi, randyall, randzder]
            loxvals_der = [randxlo, randyall, randzder]
            hiyvals_der = [randxall, randyhi, randzder]
            loyvals_der = [randxall, randylo, randzder]


            rands = [
                    hixvals_der, loxvals_der, hiyvals_der, loyvals_der,
                    ]
            randselect = random.randint(0, len(rands) - 1)
            randx = int(rands[randselect][0])
            randy = int(rands[randselect][1])
            randz = int(rands[randselect][2])

        elif layer == 3:
            # Local wound placement (interior wound region in 3D)
            randxhi = random.uniform(max(0, xboundary + lim + cell_freeze), min(nx, nx - lim - xboundary - cell_freeze))
            randxlo = random.uniform(max(0, xboundary + lim + cell_freeze), min(nx, nx - lim - xboundary - cell_freeze))
            randyall = random.uniform(max(0, yboundary + lim), min(ny, ny - lim - yboundary))
            randyhi = random.uniform(max(0, yboundary + lim + cell_freeze), min(ny, ny - lim - yboundary - cell_freeze))
            randylo = random.uniform(max(0, yboundary + lim + cell_freeze), min(ny, ny - lim - yboundary - cell_freeze))
            randxall = random.uniform(max(0, xboundary + lim), min(nx, nx - lim - xboundary))

            # Subcutaneous layer (70-100% of tissue depth)
            randzsub = random.uniform(max(0, nz - zboundary), max(1, nz - cell_freeze))


            hixvals_sub = [randxhi, randyall, randzsub]
            loxvals_sub = [randxlo, randyall, randzsub]
            hiyvals_sub = [randxall, randyhi, randzsub]
            loyvals_sub = [randxall, randylo, randzsub]

            rands = [hixvals_sub, loxvals_sub, hiyvals_sub, loyvals_sub]

            randselect = random.randint(0, len(rands) - 1)
            randx = int(rands[randselect][0])
            randy = int(rands[randselect][1])
            randz = int(rands[randselect][2])




    elif typein == 3:
        # Systemic - wound (outside wound areas in 3D)
        randxhi = random.uniform(max(0, nx - xboundary + cell_freeze), min(nx, nx - cell_freeze + lim))
        randxlo = random.uniform(max(0, cell_freeze - lim), min(nx, xboundary - cell_freeze + lim))
        randyall = random.uniform(max(0, xboundary - cell_freeze), min(ny, ny - xboundary + cell_freeze))
        randyhi = random.uniform(max(0, ny - yboundary + cell_freeze), min(ny, ny + lim - cell_freeze))
        randylo = random.uniform(max(0, cell_freeze - lim), min(ny, yboundary - cell_freeze + lim))
        randxall = random.uniform(max(0, -lim + cell_freeze), min(nx, nx + lim - cell_freeze))
        
        # Z-dimension for systemic placement (healthy tissue regions)
        # Avoid central wound layers, focus on peripheral and deep healthy tissue
        randzperipheral = random.uniform(0, max(1, zboundary - cell_freeze))  # Healthy superficial
        randzhealthy_deep = random.uniform(max(0, nz - zboundary + cell_freeze), max(1, nz - cell_freeze))  # Healthy deep
        randzavoid_wound = random.uniform(max(0, cell_freeze), max(1, zboundary - cell_freeze))  # Avoid wound core
        
        # Systemic 3D positions (avoiding wound core)
        hixvals_peri = [randxhi, randyall, randzperipheral]
        loxvals_peri = [randxlo, randyall, randzperipheral]
        hiyvals_peri = [randxall, randyhi, randzperipheral]
        loyvals_peri = [randxall, randylo, randzperipheral]
        
        hixvals_deep = [randxhi, randyall, randzhealthy_deep]
        loxvals_deep = [randxlo, randyall, randzhealthy_deep]
        hiyvals_deep = [randxall, randyhi, randzhealthy_deep]
        loyvals_deep = [randxall, randylo, randzhealthy_deep]
        
        avoid_vals = [randxall, randyall, randzavoid_wound]
        
        rands = [hixvals_peri, loxvals_peri, hiyvals_peri, loyvals_peri,
                hixvals_deep, loxvals_deep, hiyvals_deep, loyvals_deep, avoid_vals]
        
        randselect = random.randint(0, len(rands) - 1)
        randx = int(rands[randselect][0])
        randy = int(rands[randselect][1])
        randz = int(rands[randselect][2])
    
    else:
        # Default fallback - random placement anywhere in 3D space
        randx = random.randint(0, nx - 1)
        randy = random.randint(0, ny - 1)
        randz = random.randint(0, nz - 1)
    
    # Final bounds checking for all coordinates
    randx = max(0, min(int(randx), nx - 1))
    randy = max(0, min(int(randy), ny - 1))
    randz = max(0, min(int(randz), nz - 1))
    
    return [randx, randy, randz]


def safe_randopos_3D_with_types(cell_field, typein, layer, xboundary=None, yboundary=None, zboundary=None,
                                lim=None, cell_freeze=None, max_attempts=100000):
    """
    Safe 3D position finding using the type-based placement system
    
    Args:
        cell_field: CC3D cell field object
        typein: Placement type (0-3)
        Other parameters: Boundary and limit parameters
        max_attempts: Maximum placement attempts
        
    Returns:
        tuple: (x, y, z) coordinates or None if failed
    """
    
    # Get field dimensions with robust validation
    try:
        if hasattr(cell_field, 'dim'):
            dims = cell_field.dim
            nx, ny, nz = dims.x, dims.y, dims.z
        elif hasattr(cell_field, 'getDim'):
            dims = cell_field.getDim()
            nx, ny, nz = dims.x, dims.y, dims.z
        else:
            print("Could not access field dimensions - no dim attribute")
            return None
            
        # Validate dimensions
        if nx < 3 or ny < 3 or nz < 3:
            print(f"Field dimensions too small: {nx}x{ny}x{nz}")
            return None
            
    except Exception as e:
        print(f"Could not access field dimensions: {e}")
        return None
    
    # Set default parameters based on field size
    if xboundary is None:
        xboundary = nx // 4
    if yboundary is None:
        yboundary = ny // 4
    if zboundary is None:
        zboundary = nz // 4
    if lim is None:
        lim = min(10, min(nx, ny, nz) // 10)
    if cell_freeze is None:
        cell_freeze = 5
    
    # Validate parameters don't exceed field bounds
    xboundary = min(xboundary, nx // 2)
    yboundary = min(yboundary, ny // 2)
    zboundary = min(zboundary, nz // 2)
    lim = min(lim, min(nx, ny, nz) // 4)
    cell_freeze = min(cell_freeze, min(nx, ny, nz) // 8)

    # print(f"safe_randpos lim {lim}")
    # print(f"safe_randpos cell_freeze {cell_freeze}")

    for attempt in range(max_attempts):
        try:
            # Generate position using type-based strategy
            position = randopos_3D(typein, layer,xboundary, yboundary, zboundary,
                                 lim, cell_freeze, nx, ny, nz)
            
            if position is None:
                continue
                
            x, y, z = position
            
            # Validate coordinates are within bounds
            if (x < 0 or x >= nx or y < 0 or y >= ny or z < 0 or z >= nz):
                continue
            
            # Check if position is free
            if not cell_field[x, y, z]:
                return x, y, z
                
        except Exception as e:
            if attempt % 1000 == 0:
                print(f"Position search error at attempt {attempt}: {e}")
            continue
    
    # print(f"Could not find free position for typein {typein} after {max_attempts} attempts")
    return None


def safe_randopos(cell_field, typein, layer , max_attempts=1000000):
    """
    Generate a random position that is not already occupied - FIXED VERSION
    
    Args:
        cell_field: The lattice where cells are placed.
        typein: Type of position selection (passed through from original randopos).
        max_attempts: Maximum number of attempts to find an unoccupied position.
    
    Returns:
        tuple: (x, y, z) coordinates of an unoccupied position or None if failed.
    """
    
    for attempt in range(max_attempts):
        # Get position from 3D placement system
        position = safe_randopos_3D_with_types(cell_field=cell_field, typein=typein, layer=layer)
        
        # Check if position was successfully generated
        if position is None:
            # If position generation failed, try with simpler fallback
            if attempt % 10000 == 0:
                print(f"Position generation failed at attempt {attempt}, returning None")

                return None
            

        else:
            # Position was generated successfully
            randx, randy, randz = position
            
            # Double-check if position is not occupied
            if not is_position_occupied(cell_field, randx, randy, randz):
                return randx, randy, randz
    
    # If no unoccupied position is found after max_attempts
    print(f"Could not find an unoccupied position after {max_attempts} attempts")
    return None


def is_position_occupied(cell_field, x, y, z):
    """
    Check if a specific lattice position is already occupied by a cell - ROBUST VERSION
    
    Args:
        cell_field: The lattice where cells are placed.
        x (int): x-coordinate.
        y (int): y-coordinate.
        z (int): z-coordinate.
    
    Returns:
        bool: True if the position is occupied, False otherwise.
    """
    try:
        # Multiple methods to access field dimensions
        dims = None
        
        # Method 1: Standard .dim attribute
        if hasattr(cell_field, 'dim'):
            dims = cell_field.dim
        # Method 2: Alternative .getDim() method
        elif hasattr(cell_field, 'getDim'):
            dims = cell_field.getDim()
        # Method 3: Alternative .Dim attribute (capital D)
        elif hasattr(cell_field, 'Dim'):
            dims = cell_field.Dim
        else:
            print("Warning: Could not access field dimensions - using coordinate validation")
            # Try to access the position directly to test bounds
            try:
                test_cell = cell_field[x, y, z]
                # If we can access it, check if occupied
                return test_cell is not None and test_cell != 0
            except (IndexError, AttributeError):
                return True  # Treat as occupied if we can't access it
        
        # Validate coordinates are within field bounds if we have dimensions
        if dims and hasattr(dims, 'x') and hasattr(dims, 'y') and hasattr(dims, 'z'):
            if (x < 0 or x >= dims.x or y < 0 or y >= dims.y or z < 0 or z >= dims.z):
                return True  # Treat out-of-bounds as occupied
        
        # Check if the position contains a cell (not None or 0)
        cell_at_position = cell_field[x, y, z]
        return cell_at_position is not None and cell_at_position != 0
        
    except Exception as e:
        print(f"Error checking position occupation at ({x}, {y}, {z}): {e}")
        return True  # Treat error positions as occupied for safety


# Additional utility function for safe cell placement
def safe_place_cell_with_type(steppable_self, cell_type, typein, layer = 0 , cell_width=1):
    """
    Complete safe cell placement using type-based positioning
    
    Args:
        steppable_self: The steppable object containing cell_field and new_cell method
        cell_type: CC3D cell type ID
        typein: Placement type (0-3)
        cell_width: Width of cell volume
        
    Returns:
        tuple: (created_cell, x, y, z) or (None, None, None, None) if failed
    """
    
    # Get position using safe method
    position = safe_randopos(steppable_self.cell_field, typein , layer)
    
    if position is None:
        print(f"Could not place cell type {cell_type} with typein {typein}")
        return None, None, None, None
    
    x, y, z = position
    
    try:
        # Create and place cell
        new_cell = steppable_self.new_cell(cell_type)
        
        # Place cell with volume
        if cell_width == 1:
            steppable_self.cell_field[x, y, z] = new_cell
        else:
            # Validate bounds for multi-pixel placement
            dims = steppable_self.cell_field.dim
            if (x + cell_width <= dims.x and y + cell_width <= dims.y and z + cell_width <= dims.z):
                steppable_self.cell_field[x:x + cell_width, y:y + cell_width, z:z + cell_width] = new_cell
            else:
                # Use single pixel if volume doesn't fit
                steppable_self.cell_field[x, y, z] = new_cell
        
        return new_cell, x, y, z
        
    except Exception as e:
        print(f"Error placing cell: {e}")
        return None, None, None, None



def track_cell_movement(cell, dx=1.0):

    prev_x = cell.dict.get("prev_x", cell.xCOM)
    prev_y = cell.dict.get("prev_y", cell.yCOM)
    prev_z = cell.dict.get("prev_z", cell.zCOM)

    dx_step = cell.xCOM - prev_x
    dy_step = cell.yCOM - prev_y
    dz_step = cell.zCOM - prev_z

    distance = (dx_step**2 + dy_step**2 + dz_step**2)**0.5
    distance_um = distance * dx

    cell.dict["total_distance"] = cell.dict.get("total_distance", 0.0) + distance_um

    cell.dict["prev_x"] = cell.xCOM
    cell.dict["prev_y"] = cell.yCOM
    cell.dict["prev_z"] = cell.zCOM

    return distance_um

###works for 1 cell at a time
cell_movement_dict = defaultdict(list)

def track_all_cell_movements(cell_list, dx=1.0, mcs=0):


    total_movement = 0.0
    count = 0

    for cell in cell_list:
        distance_um = track_cell_movement(cell, dx)
        print(f"MCS {mcs} | Cell {cell.id} moved {distance_um:.2f}  | Total: {cell.dict['total_distance']:.2f} ") # in micrometer
        total_movement += distance_um

        cell_movement_dict[cell.id].append(distance_um)

        count += 1

    if count > 0:
        mean_movement = total_movement / count
        print(f"MCS {mcs} | Mean movement of {count} cells: {mean_movement:.2f} ") #in micrometer
    else:
        print(f" MCS {mcs} | No cells found to track.")

    print(cell_movement_dict)


def cell_death(self, cellstoDeletelist):
    cellstoDeletelist = set(cellstoDeletelist)
    cellstoDeletelist = list(cellstoDeletelist)


    deleted_cells = set()  # To track already deleted cells
    for cell in cellstoDeletelist:


        # Check if the cell is still valid (e.g., has non-zero volume)
        if cell.id in deleted_cells:
            # print(f"Cell {cell.id} has already been deleted.")
            continue  # Skip this cell as it has already been deleted

        if cell.volume == 0:
            # print(f"Cell {cell.id} has zero volume and may not exist anymore.")
            continue  # Skip this cell as it no longer exists

        # If the cell is valid and hasn't been deleted yet, proceed with deletion
        # print(f"Deleting cell {cell.id} with volume {cell.volume}")

        self.delete_cell(cell)
        print(f"the deleted cell had type {cell.type}")

        # Mark the cell as deleted
        deleted_cells.add(cell.id)

    print(r'Deleted {} cells this time'.format(len(cellstoDeletelist)))



def sigmoid(x, a, b):
    '''
    Sigmoid function:
    Returns the curve for,
    a:first inflection point
    b:second inflection point
    over x(time)
    '''
    z = np.exp(-a * (x - b))
    sig = 1 / (1 + z)

    return sig




# Get the dimensions of the lattice
nx = al.nx
ny = nx
nz = ny
dx = al.dx
dy = dx
dz = dy

L = dx * nx

# Construct the lattice in a 3D mesh
mesh = Grid3D(dx=dx, dy=dy, dz= dz, nx=nx, ny=ny, nz=nz)

# Check if the solution variable is in the mesh - necrotic cell
cellpresentnc = CellVariable(name="solution variable",
                             mesh=mesh,
                             value=0)

# Check if the solution variable is in the mesh - endothelial cell
cellpresente = CellVariable(name="solution variable",
                            mesh=mesh,
                            value=0)

# Check if the solution variable is in the mesh - platelet
cellpresentpl = CellVariable(name="solution variable",
                             mesh=mesh,
                             value=0)

# Check if the solution variable is in the mesh - keratinocyte
cellpresentkr = CellVariable(name="solution variable",
                            mesh=mesh,
                            value=0)


# Check if the solution variable is in the mesh - keratinocyte
cellpresentfib = CellVariable(name="solution variable",
                            mesh=mesh,
                            value=0)


# Check if the solution variable is in the mesh -  neutrophil
cellpresentn = CellVariable(name="solution variable",
                             mesh=mesh,
                             value=0)

# Check if the solution variable is in the mesh - apoptotic neutrophil
cellpresentap = CellVariable(name="solution variable",
                             mesh=mesh,
                             value=0)

# Check if the solution variable is in the mesh - necrotic neutrophil
cellpresentnn = CellVariable(name="solution variable",
                             mesh=mesh,
                             value=0)

# Check if the solution variable is in the mesh - mast cell
cellpresentmc = CellVariable(name="solution variable",
                             mesh=mesh,
                             value=0)

# Check if the solution variable is in the mesh - Monocytes
cellpresentmono = CellVariable(name="solution variable",
                               mesh=mesh,
                               value=0)

# Check if the solution variable is in the mesh - Macrophage 1
cellpresentm1 = CellVariable(name="solution variable",
                             mesh=mesh,
                             value=0)
# Check if the solution variable is in the mesh - Macrophage 2
cellpresentm2 = CellVariable(name="solution variable",
                             mesh=mesh,
                             value=0)

# For every solution variable disperse cytokines
cytokines = [CellVariable(name="solution variable", mesh=mesh, value=0.0) for i in range(al.total_cytokines)]




class CellSteppable(SteppableBasePy):

    def __init__(self, frequency=int(al.relaxationmcs)):
        SteppableBasePy.__init__(self, frequency)

        self.output_file = None
        self.movement_file = None
        self.movement_file_path = None
        self.cell_count_file = None

        self.summary_written = False

        # Define secretion parameters
        self.secretion_radius = 2
        self.secretion_strength = 30.0  # Increased for visibility
        self.diffusion_coeff = 1.0  # Diffusion coefficient (micrometer m²/s)
        self.decay_rate = 0.01  # Decay rate (1/s)
        self.dt = 10  # Time step for PDE solver




    def start(self):
        '''
        Initialization of simulation. Creates the field for cytokines
        and distribution of cells, with the settings.
        Also starts the plotting settings.
        '''

        global fullFileName
        global saturation_coefficient
        global endocount

        global cellpresente
        global cellpresentn


        ## creating files for tracking the data
        file_name = "simulation_data.csv"
        try:
            # The 'w' means we create a new file each time the simulation starts.
            self.output_file, full_file_path  = self.open_file_in_simulation_output_folder(file_name, "w")


            header = "mcs,ccl2_mean,ccl2_std,damps_mean,damps_std,pamps_mean,pamps_std,tgf_b_mean,tgf_b_std,pdgf_mean,pdgf_std,fgf_mean,fgf_std,tnf_a_mean,tnf_a_std,il1a_mean,il1a_std,il1b_mean,il1b_std,il6_mean,il6_std,il8_mean,il8_std,il10_mean,il10_std,il1ra_mean,il1ra_std\n"  # Add other data columns as needed
            self.output_file.write(header)
        except IOError:
            print(f"Could not open file {file_name} for writing. Aborting.")
            self.stop_simulation()

        movement_file_name = "movement_data.csv"
        try:
            self.movement_file, self.movement_file_path = self.open_file_in_simulation_output_folder(movement_file_name, "w")
            # Log the raw data: which step, which cell, what type, how far it moved.
            header = "MCS,CellID,CellType,Distance\n"
            self.movement_file.write(header)
        except IOError:
            print(f"Could not open file {movement_file_name} for writing. Aborting.")
            self.stop_simulation()

        try:
            self.cell_count_file, _ = self.open_file_in_simulation_output_folder("cell_counts.csv", "w")
            # This header is simple: MCS, the Cell's Type ID, and the Count
            header = "MCS,CellType,Count\n"
            self.cell_count_file.write(header)
        except IOError:
            print("ERROR: Could not open cell_counts.csv for writing.")




        # Create the cytokine scalar field
        self.scalarFieldil8 = self.create_scalar_field_py("il8")
        self.scalarFieldil1a = self.create_scalar_field_py("il1a")
        self.scalarFieldil1b = self.create_scalar_field_py("il1b")
        self.scalarFieldil6 = self.create_scalar_field_py("il6")
        self.scalarFieldil10 = self.create_scalar_field_py("il10")
        self.scalarFieldil1ra = self.create_scalar_field_py("il1ra")
        self.scalarFielddamps = self.create_scalar_field_py("damps")
        self.scalarFieldpamps = self.create_scalar_field_py("pamps")


        self.scalarFieldtnfa = self.create_scalar_field_py("tnfa")
        self.scalarFieldtgfb = self.create_scalar_field_py("tgfb")
        self.scalarFieldccl2 = self.create_scalar_field_py("ccl2")
        self.scalarFieldpdgf = self.create_scalar_field_py("pdgf")
        self.scalarFieldfgf = self.create_scalar_field_py("fgf")




        # placement of endothelial cells both on superficial and deeper layers (so everywhere)
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=0, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.TEMP)
            else:
                print(f"Could not place endothelial cell {i+1}/10")

        # placement of endothelial cells specifically in the dermis
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.TEMP)
            else:
                print(f"Could not place endothelial cell {i + 1}/10")

        # placement of endothelial cells specifically in the subcutaneous part
        for i in range(0, 30):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.TEMP)
            else:
                print(f"Could not place endothelial cell {i + 1}/10")

        # placement of endothelial cells specifically in the dermis outside the wound
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=3, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.TEMP)
            else:
                print(f"Could not place endothelial cell {i + 1}/10")



        # Placement of platelets on the wound edges
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=1, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.PLATELET)
            else:
                print(f"Could not place platelet cell {i+1}/10")

        # Placement of platelets within the dermis
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.PLATELET)
            else:
                print(f"Could not place platelet cell {i+1}/10")

        # Placement of platelets within the subcutaneous tissue
        for i in range(0, 30):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.PLATELET)
            else:
                print(f"Could not place platelet cell {i+1}/10")



        # Placement of keratinocytes on superficial layer
        for i in range(0, 50):
            position = safe_randopos(self.cell_field, typein=2, layer = 1)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.KERATINO)
            else:
                print(f"Could not place platelet cell {i+1}/10")


        # Necrotic cell placement both on superficial and deeper layers
        for i in range(0, 100):
            position = safe_randopos(self.cell_field, typein=0, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.NECTEMP)
            else:
                print(f"Could not place platelet cell {i+1}/10")


        # Necrotic cell placement on the wound edges
        for i in range(0, 100):
            position = safe_randopos(self.cell_field, typein=1, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.NECTEMP)
            else:
                print(f"Could not place platelet cell {i+1}/10")

        # Necrotic cell placement everywhere within the wound site
        for i in range(0, 100):
            position = safe_randopos(self.cell_field, typein=2, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.NECTEMP)
            else:
                print(f"Could not place platelet cell {i+1}/10")


        # Neutrophils within the dermis
        for i in range(0, 60):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.NEUTROPHIL)
            else:
                print(f"Could not place neutrophil {i+1}/30")


        # Neutrophil within the subcutaneous layer
        for i in range(0, 90):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.NEUTROPHIL)
            else:
                print(f"Could not place neutrophil {i+1}/30")





        # Mast cells within the dermis
        for i in range(0, 33):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.MAST)
            else:
                print(f"Could not place neutrophil {i+1}/30")

        # Mast cells within the subcutaneous layer
        for i in range(0, 32):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.MAST)
            else:
                print(f"Could not place neutrophil {i+1}/30")


        # Fibroblasts within the wound edges
        for i in range(0, 15):
            position = safe_randopos(self.cell_field, typein=1, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.FIBROBLAST)
            else:
                print(f"Could not place neutrophil {i+1}/30")


        # Fibroblasts within the dermis
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.FIBROBLAST)
            else:
                print(f"Could not place neutrophil {i+1}/30")

        # Fibroblasts within the subcutaneous layer
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.FIBROBLAST)
            else:
                print(f"Could not place neutrophil {i+1}/30")

        # monocyte on the wound border
        for i in range(0, 20):
            position = safe_randopos(self.cell_field, typein=1, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MONOCYTE)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                print(f"Could not place neutrophil {i + 1}/30")

        # monocytes within the dermis
        for i in range(0, 25):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MONOCYTE)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                print(f"Could not place neutrophil {i + 1}/30")

        # monocytes within the subcutaneous layer
        for i in range(0, 25):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MONOCYTE)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                print(f"Could not place neutrophil {i + 1}/30")


        # M1 macrophages  within the dermis
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MACROPHAGE1)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                print(f"Could not place neutrophil {i + 1}/30")


        # M1 macrophages within the subcutaneous layer
        for i in range(0, 10):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MACROPHAGE1)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                print(f"Could not place neutrophil {i + 1}/30")

        # M2 macrophages  within the dermis
        for i in range(0, 1):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MACROPHAGE2)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                print(f"Could not place neutrophil {i + 1}/30")

        # M2 macrophages within the subcutaneous layer
        for i in range(0, 1):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MACROPHAGE2)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                print(f"Could not place neutrophil {i + 1}/30")













        # All cells must have a volume and a range of movement
        for cell in self.cell_list:


            # Endothelial cell settings
            if cell.type == 1:  # Endothelial cell
                cell.dict["span"] = al.endothelial_lifespan
                cell.dict["life"] = 0
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1


            #Keratinocyte settings
            if cell.type == 10:  # Keratinocyte
                cell.dict["span"] = al.kera_lifespan
                cell.dict["life"] = random.randint(0, al.kera_lifespan)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1

            # Fibroblast settings for chemotaxis
            if cell.type == 4:  # Fibroblast
                cell.dict["span"] = al.fibro_final_lifespan
                cell.dict["life"] = al.fibro_final_lifespan
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1


            # Platelet cell settings
            if cell.type == 11:  # Platelet cell
                cell.dict["span"] = al.platelet_lifespan
                cell.dict["life"] = random.randint(0, al.platelet_lifespan)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1



            # Neutrophil settings
            if cell.type == 2:  # Neutrophil
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL8")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["life"] = random.randint(0, al.neutro_lifespan_6)


            if cell.type == 2:  # Neutrophil
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "DAMPS")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["life"] = random.randint(0, al.neutro_lifespan_6)
                # cell.dict["divide"] = vv.divpre
                # cell.dict["dividepr"] = -1

            if cell.type == 2:  # Neutrophil
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL1A")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["life"] = random.randint(0, al.neutro_lifespan_6)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1
                # cell.dict["divide"] = vv.divpre
                # cell.dict["dividepr"] = -1

            if cell.type == 2:  # Neutrophil
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL1B")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["life"] = random.randint(0, al.neutro_lifespan_6)
                # cell.dict["divide"] = vv.divpre
                # cell.dict["dividepr"] = -1


            if cell.type == 7:  # Mast Cell
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "PAMPS")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = al.mast_lifespan
                cell.dict["life"] = 0 # doesnt need specific lifespan


            # necrotic cells settings
            if cell.type == 12:
                cell.dict["span"] = al.necrotic_lifespan
                cell.dict["life"] = random.randint(0, al.necrotic_lifespan) # is it better to have a set lifespan or no?
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1

            #necrotic starter
            if cell.type == 14:
                cell.dict["span"] = al.necrotic_lifespan
                cell.dict["life"] = random.randint(0, al.necrotic_lifespan)  # is it better to have a set lifespan or no?
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1


            # endothelial starter
            if cell.type == 13:
                cell.dict["span"] = al.endothelial_lifespan
                cell.dict["life"] = random.randint(0, al.necrotic_lifespan) # is it better to have a set lifespan or no?
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1










    def step(self, mcs):
        '''
        Function defines what is done for each MCS step
        '''
        global cytokines
        global cellpresente
        global cellpresentn
        global fullFileName
        global relaxationmcs

        if mcs >= 1 and not self.summary_written:
            summary_file_name = "simulation_summary.txt"
            file_handle = None

            try:
                # Open the file and get the handle
                file_handle_tuple = self.open_file_in_simulation_output_folder(summary_file_name, "w")
                file_handle = file_handle_tuple[0]

                print(f" Writing simulation summary to {summary_file_name} ")


                file_handle.write("==== Simulation Setup ====\n")
                file_handle.write(f"Grid Size (x, y, z): {al.nx}, {al.ny}, {al.nz}\n\n")


                file_handle.write("==== Key Parameters from all_values.py ====\n")
                file_handle.write(f"Relaxation MCS: {al.relaxationmcs}\n")
                # file_handle.write(f"Calibration Factor: {getattr(vv, 'CALIBRATION_FACTOR', 'Not Defined')}\n")
                file_handle.write(f"Base Necrotic Neutrophil Chance: {al.base_necrotic_chance_neut}\n")
                file_handle.write(f"Base Necrotic Generic Cell Chance: {al.base_cell_necro_chance}\n")
                file_handle.write(f"Apoptosis Chance: {al.apop_chance_neut}\n\n")



                file_handle.write("==== Cell Counts (at MCS >= 1) ====\n")

                initial_counts = {}
                for cell in self.cell_list:
                    initial_counts[cell.type] = initial_counts.get(cell.type, 0) + 1


                type_name_map = {
                    1: "endothelial", 2: "neutrophil", 3: "monocyte",
                    4: "fibroblast", 5: "neutrophila", 6: "neutrophilnec",
                    7: "mast", 8: "macrophage1", 9: "macrophage2",
                    10: "keratino", 11: "platelet", 12: "necrotic",
                    13: "temp_endo", 14: "temp_nec"  # Use your progenitor names
                }

                for type_id, count in sorted(initial_counts.items()):
                    if count > 0:
                        type_name = type_name_map.get(type_id, f"UnknownTypeID_{type_id}")
                        file_handle.write(f"{type_name:<20}: {count}\n")

                print(" Summary file written successfully. ")

            except Exception as e:
                import traceback
                print(f"! ERROR: Could not write summary file: {e} ")
                traceback.print_exc()

            finally:
                if file_handle:
                    file_handle.close()

            self.summary_written = True






        ccount = np.zeros(al.total_celltypes + 1)
        apop_neutro_delete_list = []
        general_deletion_list = []



        # tracks the cell movement and puts the movement in a dict (for possible later calculations)
        # neutrophils = self.cell_list_by_type(self.NEUTROPHIL)
        # track_all_cell_movements(neutrophils, dx=1.0, mcs=mcs)





        # Neutrophils within the dermis
        for i in range(0, 3):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.NEUTROPHIL)
            else:
                # print(f"Could not place neutrophil {i+1}/30")
                continue


        # Neutrophil within the subcutaneous layer
        for i in range(0, 2):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.NEUTROPHIL)
            else:
                # print(f"Could not place neutrophil {i+1}/30")
                continue


        # Placement of platelets on the wound edges
        for i in range(0, 1):
            position = safe_randopos(self.cell_field, typein=1, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.PLATELET)
            else:
                # print(f"Could not place platelet cell {i+1}/10")
                continue

        # Placement of platelets within the dermis
        for i in range(0, 1):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.PLATELET)
            else:
                # print(f"Could not place platelet cell {i+1}/10")
                continue

        # Placement of platelets within the subcutaneous tissue
        for i in range(0, 2):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx+1, randy:randy+1, randz:randz+1] = self.new_cell(self.PLATELET)
            else:
                # print(f"Could not place platelet cell {i+1}/10")
                continue

        # monocyte on the wound border
        for i in range(0, 0):
            position = safe_randopos(self.cell_field, typein=1, layer=0)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MONOCYTE)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                # print(f"Could not place neutrophil {i + 1}/30")
                continue

        # monocytes within the dermis
        for i in range(0, 1):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MONOCYTE)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                # print(f"Could not place neutrophil {i + 1}/30")
                continue

        # monocytes within the subcutaneous layer
        for i in range(0, 1):
            position = safe_randopos(self.cell_field, typein=2, layer=3)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.MONOCYTE)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                # print(f"Could not place neutrophil {i + 1}/30")
                continue


        # FIBROBLAST within the dermis layer
        for i in range(0, 1):
            position = safe_randopos(self.cell_field, typein=2, layer=2)
            if position is not None:
                randx, randy, randz = position
                self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = self.new_cell(self.FIBROBLAST)
                # print("<MONOCYTES> HAVE BEEN PLACED!")
            else:
                # print(f"Could not place neutrophil {i + 1}/30")
                continue














        for cell in self.cell_list:

            # Neutrophil settings
            if cell.type == 2:  # Neutrophil
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL8")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["life"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)

            if cell.type == 2:  # Neutrophil
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "DAMPS")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["life"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                # cell.dict["divide"] = vv.divpre
                # cell.dict["dividepr"] = -1

            if cell.type == 2:  # Neutrophil
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL1A")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["life"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1
                # cell.dict["divide"] = vv.divpre
                # cell.dict["dividepr"] = -1

            if cell.type == 2:  # Neutrophil
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "IL1B")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                cell.dict["life"] = random.randint(al.neutro_lifespan_6, al.neutro_lifespan_8)
                # cell.dict["divide"] = vv.divpre
                # cell.dict["dividepr"] = -1

            # Monocyte settings for chemotaxis towards CCL2
            if cell.type == 3:  # Monocyte
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "CCL2")
                cd.setLambda(setlambda)
                cd.setChemotactTowards("MEDIUM")
                cd.setSaturationCoef(saturation_coefficient)
                cell.dict["span"] = al.mono_lifespan
                cell.dict["life"] = random.randint(0, al.mono_lifespan)


            if cell.type == 4:  # fibroblast
                cell.dict["span"] = al.fibro_final_lifespan
                cell.dict["life"] = al.fibro_final_lifespan
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1



            if cell.type == 5:  # apoptotic neutrophil
                cell.dict["span"] = al.apoptotic_neut_lifespan
                cell.dict["life"] = random.randint(0, al.apoptotic_neut_lifespan)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1


            if cell.type == 6:  # necrotic neutrophil
                cell.dict["span"] = al.necrotic_neut_lifspan
                cell.dict["life"] = random.randint(0, al.necrotic_neut_lifspan)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1

            if cell.type ==  8:  # M1 macrophage
                cell.dict["span"] = al.m1_lifespan
                cell.dict["life"] = random.randint(0, al.m1_lifespan)

            if cell.type ==  9:  # M2 macrophage
                cell.dict["span"] = al.m2_lifespan
                cell.dict["life"] = random.randint(0, al.m2_lifespan)

            if cell.type == 10:  # keratino
                cell.dict["span"] = al.kera_lifespan
                cell.dict["life"] = random.randint(0, al.kera_lifespan)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1

            if cell.type == 11:  # platelet
                cell.dict["span"] = al.platelet_lifespan
                cell.dict["life"] = random.randint(0, al.platelet_lifespan)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1

            if cell.type == 12:  # necrotic cell
                cell.dict["span"] = al.necrotic_lifespan
                cell.dict["life"] = random.randint(0, al.necrotic_lifespan)
                cell.dict["death_counter"] = cell.dict.get("death_counter", 0) + 1





        # For all cell types
        for cell in self.cell_list:
            # Save their location on the lattice
            xCOM = cell.xCOM
            yCOM = cell.yCOM
            zCOM = cell.zCOM


            xcom_before = al.nx
            xcom_after = xcom_before - 1
            # within bounds
            if xCOM >=  xcom_before:
                xCOM = xcom_after
            if yCOM >= xcom_before:
                yCOM = xcom_after
            if zCOM >= xcom_before:
                zCOM = xcom_after

            # for the positioned within bounds
            pos = int(xCOM) + (int(yCOM)) * nx + int(zCOM) * nx * ny
            # if pos > 10000:
            #     print(xCOM, yCOM)



            # endothelial
            if cell.type == 1:
                cellpresente[pos] = 1.
            # neutrophil
            if cell.type == 2:
                cellpresentn[pos] = 1.
            # monocyte
            if cell.type == 3:
                cellpresentmono[pos] = 1.
            # fibroblast
            if cell.type == 4:
                cellpresentfib[pos] = 1.
            # apoptotic neutrophil
            if cell.type == 5:
                cellpresentap[pos] = 1.
            # necrotic neutrophil
            if cell.type == 6:
                cellpresentnn[pos] = 1.
            # mast cell
            if cell.type == 7:
                cellpresentmc[pos] = 1.
            # M1 macrophage
            if cell.type == 8:
                cellpresentm1[pos] = 1.
            # M2 macrophage
            if cell.type == 9:
                cellpresentm2[pos] = 1.
            # keratinocyte
            if cell.type == 10:
                cellpresentkr[pos] = 1.
            # Platelet
            if cell.type == 11:
                cellpresentpl[pos] = 1.
            # Necrotic cells
            if cell.type == 12:
                cellpresentnc[pos] = 1.



        # Create new lists for cytokine concentrations
        ccl2_list = []
        il8_list = []
        damps_list = []
        pamps_list = []
        tgf_b_list = []
        pdgf_list = []
        fgf_list = []
        tnf_a_list = []
        il6_list = []
        il1a_list = []
        il1b_list = []
        il10_list = []
        il1ra_list = []


        tgf_beta_conc_fibro_calc = 0
        il_10_conc_fibro_calc = 0

        # For all the cell types give a position
        for cell in self.cell_list:
            xCOM = int(cell.xCOM)
            yCOM = int(cell.yCOM)
            zCOM = int(cell.zCOM)
            # zCOM = int(cell.zCOM) # SCALED removed for convenience


            xcom_before = al.nx
            xcom_after = xcom_before - 1
            # Check if cell is within boundaries
            if xCOM >= xcom_before:  # SCALED
                xCOM = xcom_after  # SCALED
            if yCOM >= xcom_before:  # SCALED
                yCOM = xcom_after
            if zCOM >= xcom_before:
                zCOM = xcom_after

                # if zCOM>=500:
                # zCOM = 499 # SCALED removed for convenience

            # Place different cytokines in the field
            ccccl2 =  self.scalarFieldccl2[xCOM, yCOM  ,zCOM]
            cil8 = self.scalarFieldil8[xCOM,yCOM,zCOM]
            cdamps = self.scalarFielddamps[xCOM,yCOM,zCOM]
            cpamps = self.scalarFieldpamps[xCOM,yCOM,zCOM]
            ctgf_b = self.scalarFieldtgfb[xCOM,yCOM,zCOM]
            cpdgf = self.scalarFieldpdgf[xCOM,yCOM,zCOM]
            cfgf = self.scalarFieldfgf[xCOM,yCOM,zCOM]
            ctnf_a = self.scalarFieldtnfa[xCOM,yCOM,zCOM]
            cil6 = self.scalarFieldil6[xCOM,yCOM,zCOM]
            cil1a = self.scalarFieldil1a[xCOM,yCOM,zCOM]
            cil1b = self.scalarFieldil1b[xCOM,yCOM,zCOM]
            cil10 = self.scalarFieldil10[xCOM,yCOM,zCOM]
            cil1ra = self.scalarFieldil1ra[xCOM,yCOM,zCOM]



            # make the values global for downward usage fibro addition chance
            tgf_beta_conc_fibro_calc = ctgf_b
            il_10_conc_fibro_calc = cil10












            
            # # Only add values that are in the tissue area
            if 10 < xCOM < 40 and 10 < yCOM < 40:

                ccl2_list.append(ccccl2)
                il8_list.append(cil8)
                damps_list.append(cdamps)
                pamps_list.append(cpamps)
                tgf_b_list.append(ctgf_b)
                pdgf_list.append(cpdgf)
                fgf_list.append(cfgf)
                tnf_a_list.append(ctnf_a)
                il6_list.append(cil6)
                il1a_list.append(cil1a)
                il1b_list.append(cil1b)
                il10_list.append(cil10)
                il1ra_list.append(cil1ra)





            # timer for death for the cells; when their life spans run out they die with a specific chance
            if cell.type == 1 or cell.type == 11 or cell.type == 4 or cell.type == 10:
                simplified_lifespan = cell.dict["span"] / 100  # because the death counter is set not at the mcs but the amount of loops a cell has lived through
                if cell.dict["death_counter"] > simplified_lifespan:
                    prob = random.random()


                    if prob <= al.base_cell_necro_chance:

                        cell.type = 12
                        cell.dict["span"] = al.necrotic_lifespan
                        cell.dict["life"] = random.randint(0, al.necrotic_lifespan)
                        cell.dict["death_counter"] = 0

                    else:
                        # print("cell became a necrotic cell")
                        # print(f"higher life then death! cell  with cell id  {cell.id} and cell type {cell.type}died at  step {mcs}")
                        general_deletion_list.append(cell)




            if cell.type == 4:
                # neutro_death_timer = cell.dict["death_counter"]
                # print(f"neutrophil with counter {neutro_death_timer}")
                prob = random.random()
                # print(prob)

                #0.01 chance to die in the first place
                if prob <= al.base_fibro_deletion_from_grid:
                    general_deletion_list.append(cell)



            ### base chance for keratinocytes to die
            if cell.type == 10:
                # neutro_death_timer = cell.dict["death_counter"]
                # print(f"neutrophil with counter {neutro_death_timer}")
                prob = random.random()
                # print(prob)

                #0.01 chance to die in the first place
                if prob <= al.base_chance_kera_death:

                    # if they die they have a chance to become necrotic or deleted where the necrotic chance is 0.01
                    # prob_2 = random.random()

                    if prob <= al.base_kera_necro_chance:
                        cell.type = 12
                        cell.dict["span"] = al.necrotic_lifespan
                        cell.dict["death_counter"] = 0
                    else:
                        general_deletion_list.append(cell)



            if cell.type == 2:


                if cell.dict["death_counter"] > cell.dict["span"] /100:

                    prob = random.random()

                    # print(prob)

                    if prob <= al.apop_chance_neut:  # prob <= 0.02
                        cell.type = 5 # Become Apoptotic Neutrophil
                        cell.dict["span"] = al.apoptotic_neut_lifespan
                        cell.dict["death_counter"] = 0
                        print(f"Neutrophil {cell.id} became APOPTOTIC at MCS {mcs}")

                    elif prob <= (al.apop_chance_neut + al.base_necrotic_chance_neut): # prob <= 0.02
                        cell.type = 6 # Become Necrotic Neutrophil
                        cell.dict["span"] = al.necrotic_neut_lifspan
                        cell.dict["death_counter"] = 0
                        print(f"Neutrophil {cell.id} became NECROTIC at MCS {mcs}")

                    else:
                        apop_neutro_delete_list.append(cell)



            # If its a resting monocyte, it transitions according to a proba    do this also for m2 but then only use Il-10 and a base proba
            if cell.type == 3:  # 7 monocyte-->8 macrophage1
                proba = (((sigmoid(ccccl2 * 10 ** 7, al.sigmoida, al.sigmoidb) * al.cccl2
                           + sigmoid(ctnf_a * 10 ** 7, al.sigmoida, al.sigmoidb) * al.lmrtnf)
                          + sigmoid(cpamps * 10 ** 7, al.sigmoida, al.sigmoidb) * al.lmrpamp)
                         - sigmoid(cil1ra * 10 ** 8, al.sigmoida, al.sigmoidb) * al.lmril1ra)

                # print(proba)
                random_value = random.random()
                if proba > random_value:
                    # print(f"Monocyte has become M1 Macrophage with probability: <{proba}> , {random_value}")
                    # Cell type changes and its characteristics are updated
                    cell.type = 8
                    cell.dict["life"] = 0
                    cell.dict["span"] = al.m1_lifespan



            # this is the main (initiation) tool for m1 to m2 macrophages
            if cell.type == 8:  # Only process M1 macrophages
                neighbor_list = self.get_cell_neighbor_data_list(cell)

                # Check if any apoptotic neutrophils (type 5) are neighboring this cell
                common_area_with_types = neighbor_list.common_surface_area_with_cell_types(cell_type_list=[5])

                if common_area_with_types != 0:
                    for neighbor, common_surface_area in neighbor_list:
                        if neighbor and neighbor.type == 5:
                            prob_m2_switch = sigmoid(cil10 * 10 ** 9, al.sigmoida, al.sigmoidb) * al.lm1il10
                            # print(f"prob m2 switch is {prob_m2_switch}")

                            # print(f"probability to switch from m1 to m2 {prob_m2_switch}")
                            if prob_m2_switch > random.random():
                                # Convert M1 to M2
                                cell.type = 9
                                cell.dict["life"] = 0
                                cell.dict["span"] = al.m2_lifespan

                                # Mark apoptotic neutrophil for deletion
                                apop_neutro_delete_list.append(neighbor)

                                print(
                                    f"M1 macrophage {cell.id} became M2 by efferocytosis of apoptotic neutrophil {neighbor.id}")
                                break  # Important: stop after first successful efferocytosis


            # this is for spontaneous switching from m1 to m2
            if cell.type == 8:  # Still an M1

                # 1. Get the local IL-10 concentration for this specific cell.
                xCOM = int(cell.xCOM)
                yCOM = int(cell.yCOM)
                zCOM = int(cell.zCOM)

                if xCOM >= al.nx: xCOM = al.nx - 1
                if yCOM >= al.ny: yCOM = al.ny - 1
                if zCOM >= al.nz: zCOM = al.nz - 1

                # (You may need your boundary checks here)
                cil10 = self.scalarFieldil10[xCOM, yCOM, zCOM]


                if cil10 > 0:
                    prob_spontaneous_switch = sigmoid(cil10, a=al.sigmoida, b=al.sigmoidb) * al.lm1il10

                    # 3. Perform the stochastic check.
                    if prob_spontaneous_switch > random.random():
                        # If the check succeeds, change the cell's type to M2.
                        cell.type = 9

                        # Optional but good practice: reset its internal age/lifespan
                        cell.dict["life"] = 0
                        cell.dict["span"] = al.m2_lifespan

                        print(
                            f"M1 macrophage {cell.id} switched to M2 spontaneously due to IL-10 (conc: {cil10:.2e}, prob: {prob_spontaneous_switch:.2f})")

            #necrotic cell removal
            if cell.type == 8 or cell.type == 9:  # M1 or M2 Macrophage


                # Check for contact with any generic necrotic cell (type 12)
                necrotic_types_to_clear = [6, 12]
                neighbor_list = self.get_cell_neighbor_data_list(cell)
                common_area_with_types = neighbor_list.common_surface_area_with_cell_types(cell_type_list=[12])

                if common_area_with_types > 0:

                    for neighbor, common_surface_area in neighbor_list:
                        # Ensure the neighbor exists and is a necrotic cell
                        if neighbor and neighbor.type in necrotic_types_to_clear:

                            # print("COMMON AREA WITH NEC CELL!")

                            # Define the probability of successful clearance

                            if al.clearance_probability_macro > random.random():


                                # Mark the necrotic cell for deletion
                                general_deletion_list.append(neighbor)

                                macrophage_type_name = "M1" if cell.type == 8 else "M2"
                                print(f"{macrophage_type_name} Macrophage {cell.id} is clearing necrotic cell {neighbor.id}")

                                # Stop this macrophage from eating more than one cell per MCS
                                break

            ccount[cell.type] += 1



        ##chance to add fibroblasts in system
        random_value = random.random()
        add_fibro_chance = sigmoid(il_10_conc_fibro_calc * 10 ** 9, al.sigmoida, al.sigmoidb) * al.lm1il10 + sigmoid(tgf_beta_conc_fibro_calc * 10 ** 9, al.sigmoida, al.sigmoidb) * al.lftgf
        if add_fibro_chance >= random_value:
            print(add_fibro_chance)

            # FIBROBLAST within the subcutaneous layer
            for i in range(0, 1):
                position = safe_randopos(self.cell_field, typein=2, layer=3)
                if position is not None:
                    randx, randy, randz = position

                    new_cell = self.new_cell(self.FIBROBLAST)
                    self.cell_field[randx:randx + 1, randy:randy + 1, randz:randz + 1] = new_cell

                    new_cell.dict["span"] = al.fibro_final_lifespan
                    new_cell.dict["life"] = al.fibro_final_lifespan
                    new_cell.dict["death_counter"] = new_cell.dict.get("death_counter", 0) + 1
                    print("<FIBRO> HAS BEEN PLACED!")
                else:
                    # print(f"Could not place neutrophil {i + 1}/30")
                    continue













        print(r'Deleting cells... mcs = {}'.format(mcs))
        # print(len(apop_neutro_delete_list))


        cell_death(self, apop_neutro_delete_list)

        cell_death(self,general_deletion_list)


        ccl2_mean = np.mean(ccl2_list)
        il8_mean = np.mean(il8_list)
        damps_mean = np.mean(damps_list)
        pamps_mean = np.mean(pamps_list)
        tgf_b_mean = np.mean(tgf_b_list)
        pdgf_mean = np.mean(pdgf_list)
        fgf_mean = np.mean(fgf_list)
        tnf_a_mean = np.mean(tnf_a_list)
        il6_mean = np.mean(il6_list)
        il1a_mean = np.mean(il1a_list)
        il1b_mean = np.mean(il1b_list)
        il10_mean = np.mean(il10_list)
        il1ra_mean = np.mean(il1ra_list)


        ccl2_std = np.std(ccl2_list)
        il8_std = np.std(il8_list)
        damps_std = np.std(damps_list)
        pamps_std = np.std(pamps_list)
        tgf_b_std = np.std(tgf_b_list)
        pdgf_std = np.std(pdgf_list)
        fgf_std = np.std(fgf_list)
        tnf_a_std = np.std(tnf_a_list)
        il6_std = np.std(il6_list)
        il1a_std = np.std(il1a_list)
        il1b_std = np.std(il1b_list)
        il10_std = np.std(il10_list)
        il1ra_std = np.std(il1ra_list)

        #data writing fro the movement the count of cells and the amount of cytokine secreted

        if self.movement_file:
            # Iterate through every single cell in the simulation
            for cell in self.cell_list:
                # Use your excellent function to get the distance moved in this step
                distance_this_step = track_cell_movement(cell, dx=al.dx)

                # Get other important info
                cell_id = cell.id
                cell_type = cell.type

                # Format the data into a CSV line
                data_string = f"{int(mcs)},{int(cell_id)},{int(cell_type)},{float(distance_this_step)}\n"

                # Write it to the file
                self.movement_file.write(data_string)

        if self.cell_count_file:
            # Create a dictionary to store counts, initialized to zero
            counts = {i: 0 for i in range(1, al.total_celltypes + 20)}  # a bit more then the total cell types
            # Loop through cells once to get counts
            for cell in self.cell_list:
                if cell.type in counts:
                    counts[cell.type] += 1

            # Write one line for each cell type that has a count > 0
            for cell_type, count in counts.items():
                if count > 0:
                    data_string = f"{mcs},{cell_type},{count}\n"
                    self.cell_count_file.write(data_string)


        data_line_to_write = f"{mcs},{ccl2_mean},{ccl2_std},{damps_mean},{damps_std},{pamps_mean},{pamps_std},{tgf_b_mean},{tgf_b_std},{pdgf_mean},{pdgf_std},{fgf_mean},{fgf_std},{tnf_a_mean},{tnf_a_std},{il1a_mean},{il1a_std},{il1b_mean},{il1b_std},{il6_mean},{il6_std},{il8_mean},{il8_std},{il10_mean},{il10_std},{il1ra_mean},{il1ra_std}\n"
        if self.output_file:
            self.output_file.write(data_line_to_write)
            # Optional but good practice: flush the buffer to make sure it's written to disk
            self.output_file.flush()



        print(f"The means of the following cytokines are: DAMPS {damps_mean}, PAMPS {pamps_mean}, IL8 {il8_mean}, CCL2 {ccl2_mean}")
        print(f"The mean of TGF-B is {tgf_b_mean}, for PDGF is {pdgf_mean} for FGF is {fgf_mean}")
        print(f"The mean of TNF-A is {tnf_a_mean}, for IL-6 is {il6_mean}, for IL-1A is {il1a_mean}, for IL-1B is {il1b_mean}")
        print(f"The mean of IL-10 is {il10_mean}, for IL-1RA is {il1ra_mean}")






        if self.output_file:
            self.output_file.flush()
        if self.movement_file:
            self.movement_file.flush()
        if self.cell_count_file:
            self.cell_count_file.flush()
    def finish(self):
        """
        Finish Function is called after the last MCS
        """
        if self.output_file:
            self.output_file.close()

        if self.movement_file:
            self.movement_file.close()

        if self.cell_count_file:
            self.cell_count_file.close()



    def on_stop(self):
        # this gets called each time user stops simulation

        if self.output_file:
            self.output_file.close()


        if self.movement_file:
            self.movement_file.close()
        # return

        if self.cell_count_file:
            self.cell_count_file.close()



class CellDifferentiationSteppable(SteppableBasePy):
    def __init__(self, frequency=1, transition_mcs=10):
        SteppableBasePy.__init__(self, frequency)
        self.transition_mcs = transition_mcs
        self.have_transitioned = False

    def step(self, mcs):
        # We only want this to run once, at the specified MCS.
        if mcs == self.transition_mcs and not self.have_transitioned:
            print(f" Running Differentiation at MCS={mcs}")

            # This loop finds every 'movable_cell' and differentiates it into 'endothelial'.
            end_cell_count = 0
            nec_cell_count = 0
            for cell in self.cell_list:

                if cell.type == 13:
                    cell.type = 1
                    end_cell_count += 1
                    # print(f"     - Differentiated {end_cell_count} movable_cells into endothelial cells.")

                if cell.type == 14:
                    cell.type = 12
                    nec_cell_count += 1
                    # print(f"     - Differentiated {nec_cell_count} movable_cells into necrotic cells.")



            self.have_transitioned = True






