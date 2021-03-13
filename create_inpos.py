def longinitp(directory, numdays):    

    txt = ''    

    for i in range( 1, 24*(numdays) + 1):
        editp = {'p': i}
        txt = '\n'.join( 
                (txt, 
                "305 7.5 24.9 {p} 1.0\n305 350 26.4 {p} 1.0\n305 350 27.7 {p} 1.0".format( **editp ))
                )

    file = open (directory + "/initial_positions.txt", "w")
    file.write(txt)
    
    print ("initial_positions.txt created in ", directory)
    
    # the indices for depth correspond to the w grid (in variable gdepw_0 in meshmask file)
        
    # 303, 445, 24.9 = 50m
    # 303, 445, 26.4 = 75m
    # 303, 445, 27.7 = 100m
    # 305, 350, * = Boundary Pass(ish)
    # 124, 737, 24.9 = Discovery Pass @ 50m
