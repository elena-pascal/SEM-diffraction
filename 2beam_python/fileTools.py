

def readInput(fileName):
    '''
    Read parameters from a given file.
    '''
    with open(fileName) as file:
        #split into lines
        rawdata = file.read().split('\n')

    data = {}

    for line in rawdata:
        if line.strip(): # ignore empty lines
            if not line.startswith("#"): # ignore comment lines
                parameter = line.split(":")

                # trim whitespaces
                param = [p.strip() for p in parameter]
                # g vector is a list
                if param[0] == 'g':
                    gv = param[1].split(',')
                    data['g'] = [float(g.strip()) for g in gv]

                elif param[0] == 'type':
                    data[param[0]] = param[1]

                else:
                    # assign to dictionary
                    data[param[0]] = float(param[1])

    return data

def fileName(intData, folder):

    parameters = indata['type'] + '-' + str(int(indata['tiltS'])) + '_' + \
                    str(abs(indata['tiltS'] - int(indata['tiltS'])))[2:] + 'tilt' +  \
                    '-' + str(copysign(int(wi),wi))[0:2] + \
                    '_' + str(abs(wi - int(wi)))[2:] + 'w'

    bright_field = folder + '/perfect_bright' + '-' + parameters  + '.out'
    dark_field = folder + '/perfect_dark' + '-' + parameters  + '.out'
    contrast_file = folder + '/perfect_contrast' + '-' + parameters  + '.out'
    depth_file = folder + '/perfect_depth' + '-' + parameters + '.out'

    return [bright_field, dark_field, contrast_file, depth_file]



def writeOutput(contrast, fileNameC, depth, fileNameD):
    '''
    Write the contrast and depth map data to separate files
    '''
    with open(fileNameC, 'w+') as file:
        for row in contrast:
            for item in row:
                file.write('%s\n' % item)

    with open(fileNameD,'w+') as file:
        for row in depth:
            for item in row:
                file.write('%s\n' % item)
