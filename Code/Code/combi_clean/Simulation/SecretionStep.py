from cc3d.core.PySteppables import SteppableBasePy
import numpy as np
import all_values as al
from diffusion_calculation import diffusion_calculation

import __main__ as main

nx = al.nx
ny = nx
nz = ny
dx = al.dx
dy = dx
dz = dy

L = dx * nx



class SecretionSteppables(SteppableBasePy):

    def __init__(self, frequency=500):
        SteppableBasePy.__init__(self, frequency)

        # Define secretion parameters
        self.secretion_radius = 2
        self.secretion_strength = 1.0  # Increased for visibility
        self.diffusion_coeff = 1.0  # Diffusion coefficient
        self.decay_rate = 0.01  # Decay rate (1/s)
        self.dt = 10  # Time step for PDE solver

        # Pre-calculate neighbor offsets
        self.neighbor_offsets = []
        for dx_off in range(-self.secretion_radius, self.secretion_radius + 1):
            for dy_off in range(-self.secretion_radius, self.secretion_radius + 1):
                for dz_off in range(-self.secretion_radius, self.secretion_radius + 1):
                    distance = np.sqrt(dx_off ** 2 + dy_off ** 2 + dz_off ** 2)
                    if distance <= self.secretion_radius:
                        self.neighbor_offsets.append((dx_off, dy_off, dz_off, distance))

        print(f"THESE ARE THE NEIGHBOR OFFSETS!{self.neighbor_offsets}")

    def cytokine_source_term(self):

        # This is more direct than creating a 3D array and flattening it.
        source_term_CCL2 = np.zeros(al.nx * al.ny * al.nz)
        source_term_IL8 = np.zeros(al.nx * al.ny * al.nz)
        source_term_DAMPS = np.zeros(al.nx * al.ny * al.nz)
        source_term_PAMPS = np.zeros(al.nx * al.ny * al.nz)

        source_term_TGF_B = np.zeros(al.nx * al.ny * al.nz)
        source_term_PDGF = np.zeros(al.nx * al.ny * al.nz)
        source_term_FGF = np.zeros(al.nx * al.ny * al.nz)
        source_term_TNF_A = np.zeros(al.nx * al.ny * al.nz)

        source_term_IL6 = np.zeros(al.nx * al.ny * al.nz)
        source_term_IL1A = np.zeros(al.nx * al.ny * al.nz)
        source_term_IL1B = np.zeros(al.nx * al.ny * al.nz)

        source_term_IL10 = np.zeros(al.nx * al.ny * al.nz)
        source_term_IL1RA = np.zeros(al.nx * al.ny * al.nz)

        for cell in self.cell_list:
            xCOM = int(np.clip(cell.xCOM, 0, nx - 1))
            yCOM = int(np.clip(cell.yCOM, 0, ny - 1))
            zCOM = int(np.clip(cell.zCOM, 0, nz - 1))

            # Add secretion in neighborhood around cell
            for dx_off, dy_off, dz_off, distance in self.neighbor_offsets:
                x = xCOM + dx_off
                y = yCOM + dy_off
                z = zCOM + dz_off
                # Check boundaries

                if 0 <= x < nx and 0 <= y < ny and 0 <= z < nz:
                    # Convert 3D coordinates to mesh cell ID
                    cell_id = x + y * nx + z * nx * ny



                    if cell.type == self.ENDOTHELIAL:
                        secretion_amount_CCL2 = self.secretion_strength * (1 + distance * 0.01)
                        source_term_CCL2[cell_id] += secretion_amount_CCL2

                    elif cell.type == self.PLATELET:
                        secretion_amount_TGF_B = self.secretion_strength * (1 + distance * 0.01)
                        source_term_TGF_B[cell_id] += secretion_amount_TGF_B

                        secretion_amount_FGF = self.secretion_strength * (1 + distance * 0.01)
                        source_term_FGF[cell_id] += secretion_amount_FGF

                        secretion_amount_PDGF = self.secretion_strength * (1 + distance * 0.01)
                        source_term_PDGF[cell_id] += secretion_amount_PDGF


                    elif cell.type == self.KERATINO:
                        secretion_amount_CCL2 = self.secretion_strength *(1 + distance * 0.01)
                        source_term_CCL2[cell_id] += secretion_amount_CCL2

                        # if secretion_amount_CCL2 > 0:

                    elif cell.type == self.FIBROBLAST:
                        secretion_amount_CCL2 = self.secretion_strength * (1 + distance * 0.01)
                        source_term_CCL2[cell_id] += secretion_amount_CCL2

                    elif cell.type == self.MAST:
                        secretion_amount_IL8 = self.secretion_strength * (1 + distance * 0.01)
                        source_term_IL8[cell_id] += secretion_amount_IL8

                        secretion_amount_IL6 = self.secretion_strength * (1 + distance * 0.01)
                        source_term_IL6[cell_id] += secretion_amount_IL6

                        secretion_amount_TNF_A = self.secretion_strength * (1 + distance * 0.01)
                        source_term_TNF_A[cell_id] = secretion_amount_TNF_A


                    elif cell.type == self.NECROTIC:
                        secretion_amount_DAMPS = self.secretion_strength * (1 + distance * 0.01)
                        source_term_DAMPS[cell_id] += secretion_amount_DAMPS

                        secretion_amount_PAMPS = self.secretion_strength * (1 + distance * 0.01)
                        source_term_PAMPS += secretion_amount_PAMPS


                    elif cell.type == self.NEUTROPHIL:
                        secretion_amount_IL8 = self.secretion_strength * (1 + distance * 0.01)
                        source_term_IL8[cell_id] += secretion_amount_IL8

                        secretion_amount_IL6 = self.secretion_strength * (1 + distance * 0.01)
                        source_term_IL6[cell_id] += secretion_amount_IL6

                        secretion_amount_TNF_A = self.secretion_strength * (1 + distance * 0.01)
                        source_term_TNF_A[cell_id] = secretion_amount_TNF_A

                        secretion_amount_IL1B = self.secretion_strength * (1 + distance * 0.01)
                        source_term_IL1B += secretion_amount_IL1B


                    elif cell.type == self.NEUTROPHILNEC:
                        secretion_amount_DAMPS = self.secretion_strength * (1 + distance * 0.01)
                        source_term_DAMPS[cell_id] += secretion_amount_DAMPS

                    elif cell.type == self.MACROPHAGE1:
                        secretion_amount_IL1A = self.secretion_strength * (1 + distance * 0.01)
                        source_term_IL1A += secretion_amount_IL1A

                        secretion_amount_IL1B = self.secretion_strength * (1 + distance * 0.01)
                        source_term_IL1B += secretion_amount_IL1B

                    elif cell.type == self.MACROPHAGE2:
                        secretion_amount_TGF_B = self.secretion_strength * (1 + distance * 0.01)
                        source_term_TGF_B += secretion_amount_TGF_B

                        secretion_amount_PDGF = self.secretion_strength* (1 + distance * 0.01)
                        source_term_PDGF += secretion_amount_PDGF

                        secretion_amount_FGF = self.secretion_strength * (1 + distance * 0.01)
                        source_term_FGF += secretion_amount_FGF

                        secretion_amount_IL10 = self.secretion_strength* (1 + distance * 0.01)
                        source_term_IL10 += secretion_amount_IL10

                        secretion_amount_Il1RA = self.secretion_strength * (1 + distance * 0.01)
                        source_term_IL1RA += secretion_amount_Il1RA

        return [source_term_CCL2, source_term_IL8, source_term_DAMPS, source_term_PAMPS, source_term_TGF_B,
                source_term_PDGF, source_term_FGF, source_term_TNF_A, source_term_IL6, source_term_IL1A,
                source_term_IL1B, source_term_IL10, source_term_IL1RA]






    def start(self):


        self.scalarFieldil8 = self.field.il8
        self.scalarFieldil1a = self.field.il1a
        self.scalarFieldil1b = self.field.il1b
        self.scalarFieldil6 = self.field.il6
        self.scalarFieldil10 = self.field.il10
        self.scalarFieldil1ra = self.field.il1ra
        self.scalarFielddamps = self.field.damps
        self.scalarFieldpamps = self.field.pamps


        self.scalarFieldtnfa = self.field.tnfa
        self.scalarFieldtgfb = self.field.tgfb
        self.scalarFieldccl2 = self.field.ccl2
        self.scalarFieldpdgf = self.field.pdgf
        self.scalarFieldfgf = self.field.fgf


    def step(self, mcs):
        # Called every `frequency` MCS
        print(f"MyNewSteppable running at MCS={mcs}")
        # Example action: increment some global array or cytokine


        source_term_CCL2, source_term_IL8, source_term_DAMPS, source_term_PAMPS, source_term_TGF_B, source_term_PDGF, source_term_FGF, source_term_TNF_A, source_term_IL6, source_term_IL1A, source_term_IL1B, source_term_IL10, source_term_IL1RA = self.cytokine_source_term()
        main.cytokines = diffusion_calculation(main.cellpresente, main.cellpresentn, main.cellpresentnc, main.cellpresentpl, main.cellpresentkr, main.cellpresentfib,
                                               main.cellpresentmc, main.cellpresentnn, main.cellpresentm1, main.cellpresentm2, source_term_CCL2,
                                               source_term_IL8, source_term_DAMPS, source_term_PAMPS, source_term_TGF_B, source_term_PDGF,
                                               source_term_FGF, source_term_TNF_A, source_term_IL6, source_term_IL1A, source_term_IL1B,
                                               source_term_IL10, source_term_IL1RA, main.cytokines, main.mesh)

        # Create the scalr field to place the cytokines

        self.scalarFieldccl2[:] = np.reshape(main.cytokines[0], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldil8[:] = np.reshape(main.cytokines[1], (al.nx, al.ny, al.nz), 'F')
        self.scalarFielddamps[:] = np.reshape(main.cytokines[2], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldpamps[:] = np.reshape(main.cytokines[3], (al.nx, al.ny, al.nz), 'F')

        self.scalarFieldtgfb[:] = np.reshape(main.cytokines[4], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldpdgf[:] = np.reshape(main.cytokines[5], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldfgf[:] = np.reshape(main.cytokines[6], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldtnfa[:] = np.reshape(main.cytokines[7], (nx, ny, nz), 'F')

        self.scalarFieldil6[:] = np.reshape(main.cytokines[8], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldil1a[:] = np.reshape(main.cytokines[9], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldil1b[:] = np.reshape(main.cytokines[10], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldil10[:] = np.reshape(main.cytokines[11], (al.nx, al.ny, al.nz), 'F')
        self.scalarFieldil1ra[:] = np.reshape(main.cytokines[12], (al.nx, al.ny, al.nz), 'F')
