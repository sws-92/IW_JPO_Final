def namelist(directory, n, tl, lmt = None):
    
    indn = n
    
    input_digits = len(str(n))
    
    ind0 = "1"
    
    while len(ind0) < input_digits:
        ind0 = "0" + ind0
        
    
        
        

        
    edit = {"dir": directory, "lmt": n*24, "output": 3*24*tl-1, "input_digits": input_digits, "ind0": ind0, "indn": indn, "maxsize": 3 }
    if lmt != None:
        edit = {"dir": directory, "length": n, "lmt": lmt, "output": tl*24*3, "input_digits": input_digits, "ind0": ind0, "indn": indn, "maxsize": 3}
    txt = """&ARIANE
        key_alltracers =.FALSE.,
        key_sequential =.TRUE.,
        key_ascii_outputs =.TRUE.,
        mode ='qualitative',
        forback ='forward',
        bin ='nobin',
        init_final ='init',
        nmax = 999999,
        tunit = 3600,
        ntfic = 1,
        tcyc =86400,

/

&OPAPARAM
        imt =398,
        jmt =898,
        kmt =40,
        lmt = {lmt},
        key_periodic =.FALSE.,
        key_jfold =.FALSE.,
        key_computew =.FALSE.,
        key_partialsteps =.TRUE.,
/

&QUALITATIVE
        delta_t = 1200,
        frequency = 1,
        nb_output = {output}, 
        key_region =.FALSE.,
/

&ZONALCRT
        c_dir_zo ='{dir}',
        c_prefix_zo ='SalishSea_',
        ind0_zo ={ind0},
        indn_zo ={indn},
        maxsize_zo ={maxsize},
        c_suffix_zo ='_grid_U.nc',
        nc_var_zo ='vozocrtx',
        nc_var_eivu ='NONE',
        nc_att_mask_zo ='NONE',
/

&MERIDCRT
        c_dir_me ='{dir}',
        c_prefix_me ='SalishSea_',
        ind0_me ={ind0},
        indn_me ={indn},
        maxsize_me ={maxsize},
        c_suffix_me ='_grid_V.nc',
        nc_var_me ='vomecrty',
        nc_var_eivv ='NONE',
        nc_att_mask_me ='NONE',
/

&VERTICRT
	c_dir_ve     = '{dir}',
	c_prefix_ve  = 'SalishSea_',
	ind0_ve      = {ind0},
	indn_ve      = {indn},
	maxsize_ve   = {maxsize},
	c_suffix_ve  = '_grid_W.nc',
	nc_var_ve    = 'vovecrtz',
	nc_var_eivw  = 'NONE',
	nc_att_mask_ve = 'NONE',
/

&MESH
        dir_mesh ='/ocean/eolson/MEOPAR/NEMO-forcing/grid',
        fn_mesh ='mesh_mask201702.nc',
        nc_var_xx_tt ='glamt',
        nc_var_xx_uu ='glamu',
        nc_var_yy_tt ='gphit',
        nc_var_yy_vv ='gphiv',
        nc_var_zz_ww ='gdepw_0',
        nc_var_e2u ='e2u',
        nc_var_e1v ='e1v',
        nc_var_e1t ='e1t',
        nc_var_e2t ='e2t',
        nc_var_e3t ='e3t_0',
        nc_var_tmask ='tmask',
        nc_mask_val =0.,
/


&SEQUENTIAL
	maxcycles =1,
/
""".format(**edit)
    
    with open(directory + "/namelist", 'w') as file:
        file.write(txt)


