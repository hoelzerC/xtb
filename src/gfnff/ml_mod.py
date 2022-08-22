def receive_ml_input_send_output(xyz):
    import numpy as np
    ###################################
    #   replace this with ML stuff    #
    ###################################
    # define nonsense energy term
    energy = np.sum(abs(xyz)) # just sum up coordinates
    print('_python: Energy=%25.15f' % energy)
    # define nonsense gradient
    gradient = xyz*0.1 # 0.1 times atom position
    ###################################
    #       replace up to here        #
    ###################################
    # define dictionary containing the energy and gradient
    return_dict = {}
    return_dict['energy']=energy
    return_dict['gradient']=gradient
    # return energy and gradient in dictionary
    return return_dict

