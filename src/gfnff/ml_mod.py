def receive_ml_input_send_output(xyz):
    # define nonsense energy term
    energy = float(sum(sum(xyz, []))) # just sum up coordinates
    # define nonsense gradient
    gradient = float(xyz*0.1) # 0.1 times atom position

    #return (energy, gradient) # consider using dict, list or class here
    return (energy) # first only return energy

def testFct():
    # just plain print
    print('_python: This is testFct().')
