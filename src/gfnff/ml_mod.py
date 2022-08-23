def receive_ml_input_send_output(*args, **kwargs):
    import numpy as np
    ###################################
    #   replace this with ML stuff    #
    ###################################
    # define nonsense energy term
    xyz=np.transpose(kwargs["xyz"])
    eatoms=kwargs["eatoms"]
    energy = np.sum(abs(xyz)) # just sum up coordinates
    gradient=xyz
    for i in range(3):
        for j in range(len(eatoms)):
            gradient[j,i] = xyz[j,i]*eatoms[j] # 0.1 times atom position
    ###################################
    #       replace up to here        #
    ###################################
    # define dictionary containing the energy and gradient
    return_dict = {}
    return_dict['energy']=energy
    return_dict['gradient']=gradient
    # return energy and gradient in dictionary
    return return_dict

def testFct():
    # just plain print
    print('_python: This is testFct().')

def testPythonsImport():
    import os
    #import numpy
    print("PYTHONPATH:", os.environ.get('PYTHONPATH'))
    #print("PATH:", os.environ.get('PATH'))
    #print('Imported numpy version:', numpy.version.version)
