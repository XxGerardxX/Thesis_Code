from fipy import CellVariable, TransientTerm, DiffusionTerm, ImplicitSourceTerm, LinearGMRESSolver

from all_values import *

def diffusion_calculation(cellpresente, cellpresentn, cellpresentnc, cellpresentpl, cellpresentkr, cellpresentfib, cellpresentmc, cellpresentnn, cellpresentm1, cellpresentm2, source_term_CCL2, source_term_IL8, source_term_DAMPS, source_term_PAMPS, source_term_TGF_B, source_term_PDGF, source_term_FGF, source_term_TNF_A, source_term_IL6, source_term_IL1A, source_term_IL1B, source_term_IL10, source_term_IL1RA, cytokines, mesh):
    '''
    Solves the diffusion equations with a mesh into consideration
    '''



    # Define the Variables being solved
    ccl2 = cytokines[0]
    il8 = cytokines[1]
    damps = cytokines[2]
    pamps = cytokines[3]

    tgf_b = cytokines[4]
    pdgf = cytokines[5]
    fgf = cytokines[6]
    tnf_a = cytokines[7]
    il6 = cytokines[8]
    il1a = cytokines[9]
    il1b = cytokines[10]
    il10 = cytokines[11]
    il1ra = cytokines[12]





    ccl2 = CellVariable(mesh = mesh,value= ccl2)
    il8 = CellVariable(mesh = mesh,value= il8)
    damps = CellVariable(mesh = mesh,value= damps)
    pamps = CellVariable(mesh = mesh,value= pamps)

    tgf_b = CellVariable(mesh = mesh,value= tgf_b)
    pdgf = CellVariable(mesh = mesh,value= pdgf)
    fgf = CellVariable(mesh = mesh,value= fgf)
    tnf_a = CellVariable(mesh = mesh,value= tnf_a)

    il6 = CellVariable(mesh = mesh,value = il6)
    il1a = CellVariable(mesh = mesh,value = il1a)
    il1b = CellVariable(mesh=mesh, value=il1b)
    il10 = CellVariable(mesh=mesh, value=il10)
    il1ra = CellVariable(mesh=mesh, value=il1ra)



    # Test with explicit Dirichlet boundaries on all faces
    ccl2.constrain(0., mesh.exteriorFaces)
    il8.constrain(0.,mesh.exteriorFaces)
    damps.constrain(0., mesh.exteriorFaces)
    pamps.constrain(0., mesh.exteriorFaces)

    tgf_b.constrain(0., mesh.exteriorFaces)
    pdgf.constrain(0., mesh.exteriorFaces)
    fgf.constrain(0., mesh.exteriorFaces)
    tnf_a.constrain(0., mesh.exteriorFaces)

    il6.constrain(0., mesh.exteriorFaces)
    il1a.constrain(0., mesh.exteriorFaces)
    il1b.constrain(0., mesh.exteriorFaces)
    il10.constrain(0., mesh.exteriorFaces)

    il1ra.constrain(0., mesh.exteriorFaces)

    # create extra variable to add source_term
    source_var = CellVariable(mesh=mesh, value=source_term_CCL2)
    source_var_IL8 = CellVariable(mesh=mesh, value=source_term_IL8)
    source_var_DAMPS = CellVariable(mesh=mesh, value=source_term_DAMPS)
    source_var_PAMPS = CellVariable(mesh=mesh, value=source_term_PAMPS)

    source_var_TGF_B =CellVariable(mesh=mesh, value=source_term_TGF_B)
    source_var_PDGF =CellVariable(mesh=mesh, value=source_term_PDGF)
    source_var_FGF = CellVariable(mesh=mesh, value=source_term_FGF)
    source_var_TNF_A = CellVariable(mesh=mesh, value=source_term_TNF_A)

    source_var_IL6 = CellVariable(mesh=mesh, value=source_term_IL6)
    source_var_IL1A = CellVariable(mesh=mesh, value=source_term_IL1A)
    source_var_IL1B = CellVariable(mesh=mesh, value=source_term_IL1B)
    source_var_IL10 = CellVariable(mesh=mesh, value=source_term_IL10)

    source_var_IL1RA = CellVariable(mesh=mesh, value=source_term_IL1RA)






    mysolver=LinearGMRESSolver()






    ## ccl2 test equation using il8
    #TODO FIX THE IMPLICIT SOURCE TERMPS AS WELL AND EXTEND THE FORMULAS

    eqccl2= TransientTerm() == DiffusionTerm(coeff=D_ccl2) - ImplicitSourceTerm(mu_ccl2) + (k_endoccl2 * cellpresente) * source_var  + (k_kera_ccl2 * cellpresentkr) * source_var + (k_fibro_ccl2 * cellpresentfib) * source_var
    eqil8 = TransientTerm() == DiffusionTerm(coeff=D_il8) - ImplicitSourceTerm(mu_il8) + (k_neutroil8 * cellpresentn) * source_var_IL8 + (k_mastil8 * cellpresentmc) * source_var_IL8
    eqdamps = TransientTerm() == DiffusionTerm(coeff=D_DAMPS) - ImplicitSourceTerm(mu_DAMPS)   + (k_neutro_nec_DAMPS * cellpresentnn) * source_var_DAMPS + (k_nec_DAMPS * cellpresentnc) * source_var_DAMPS


    eqpamps = TransientTerm() == DiffusionTerm(coeff=D_PAMPS) - ImplicitSourceTerm(mu_PAMPS) + (k_nec_PAMPS * cellpresentnc) * source_var_PAMPS

    eqtgf_b = TransientTerm() == DiffusionTerm(coeff=D_tgfb) - ImplicitSourceTerm(mu_tgfb) + (k_m2_tgfb * cellpresentm2) * source_var_TGF_B + (k_plat_tgfb * cellpresentpl)  * source_var_TGF_B
    eqpdgf = TransientTerm() == DiffusionTerm(coeff=D_PDGF) - ImplicitSourceTerm(mu_PDGF) + (k_plat_PDGF * cellpresentpl) * source_var_PDGF + (k_m2_PDGF * cellpresentm2) * source_var_PDGF
    eqfgf = TransientTerm() == DiffusionTerm(coeff=D_FGF) - ImplicitSourceTerm(mu_FGF) + (k_plat_FGF * cellpresentpl) * source_var_FGF + (k_m2_FGF * cellpresentm2) * source_var_FGF
    eqtnf_a = TransientTerm() == DiffusionTerm(coeff=D_tnfa) - ImplicitSourceTerm(mu_tnfa) + (k_mast_tnfa * cellpresentmc) * source_var_TNF_A + (k_neut_tnfa * cellpresentn) * source_var_TNF_A

    eqil6 = TransientTerm() == DiffusionTerm(coeff=D_il6) - ImplicitSourceTerm(mu_il6) + (k_mast_il6 * cellpresentmc) * source_var_IL6 + (k_neut_il6 * cellpresentn) * source_var_IL6
    eq_il1a = TransientTerm() == DiffusionTerm(coeff=D_il1a) - ImplicitSourceTerm(mu_il1a) + (k_m1_il1a * cellpresentm1) * source_var_IL1A
    eq_il1b = TransientTerm() == DiffusionTerm(coeff=D_il1b) - ImplicitSourceTerm(mu_il1b) + (k_m1_il1b * cellpresentm1) * source_var_IL1B  + (k_neut_il1b * cellpresentn) * source_var_IL1B
    eqil10 = TransientTerm() == DiffusionTerm(coeff=D_il10) - ImplicitSourceTerm(mu_il10) + (k_m2_il10 * cellpresentm2) * source_var_IL10

    eqil1ra = TransientTerm() == DiffusionTerm(coeff=D_il1RA) - ImplicitSourceTerm(mu_il1RA) + (k_m2_il1RA * cellpresentm2) * source_var_IL1RA




    # Solve for the duration(dt)
    for i in range(fipy_duration):


        ## ccl2 test solve using il8
        eqccl2.solve(var=ccl2, dt=1.0, solver=mysolver)

        ##il8 test solve
        eqil8.solve(var=il8, dt=1.0, solver=mysolver)

        ## damps test
        eqdamps.solve(var=damps, dt=1.0, solver=mysolver)

        ## pamps test
        eqpamps.solve(var=pamps, dt=1.0, solver=mysolver)

        ##tgf_b test
        eqtgf_b.solve(var=tgf_b, dt=1.0, solver=mysolver)

        ## pdgf test
        eqpdgf.solve(var=pdgf, dt=1.0, solver=mysolver)

        ## fgf test
        eqfgf.solve(var=fgf, dt=1.0, solver=mysolver)

        #tnf_a test
        eqtnf_a.solve(var=tnf_a, dt=1.0, solver=mysolver)

        #il-6 test
        eqil6.solve(var=il6, dt=1.0, solver=mysolver)

        #il1a test
        eq_il1a.solve(var=il1a, dt=1.0, solver=mysolver)

        #il1b test
        eq_il1b.solve(var=il1b, dt=1.0, solver=mysolver)

        #il10 test
        eqil10.solve(var=il10, dt=1.0, solver=mysolver)

        #il1ra solve
        eqil1ra.solve(var=il1ra, dt=1.0, solver=mysolver)



    return [ccl2,il8,damps,pamps, tgf_b, pdgf, fgf, tnf_a, il6, il1a, il1b, il10, il1ra]



