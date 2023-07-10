from pathlib import Path

from scipy.misc import electrocardiogram

area = ['$16^\smblksquare$', '$25^\smblksquare$', '$35^\smblksquare$',
        '$50^\smblksquare$', '$70^\smblksquare$', '$95^\smblksquare$', '$120^\smblksquare$', '$150^\smblksquare$', '$185^\smblksquare$']


def plotingLines(data, location, name, lineSize):
    dir_path = Path(location)
    file_name = name
    nl = data['A'].shape[1]
    lineSize = lineSize.astype(int)
    k = open(dir_path.joinpath(file_name), 'w')
    k.write('u, v, lw, style, color, label\n')
    for i in range(0, nl):
        f = data['f'][i]+1
        t = data['t'][i]+1
        k.write('%.0d, %.0d, 1, solid, black, %s\n' %
                (f, t, area[lineSize[i]]))
    k.close()


def plotingBuses(data, location, name):
    dir_path = Path(location)
    file_name = name
    nb = data['A'].shape[0]
    xx, yy = data['xy'].T
    k = open(dir_path.joinpath(file_name), 'w')
    k.write('id, x, y, size, style, color, label\n')
    for i in range(0, nb):
        id = i+1
        x = xx[i]
        y = yy[i]
        label = i+1
        if id == 1:
            k.write('%.0f, %.2f, %.2f, 0.3, thin, gray, %.0f\n' %
                    (id, x, y, label))
        else:
            k.write('%.0f, %.2f, %.2f, 0.25, thin, white, %.0f\n' %
                    (id, x, y, label))
    k.close()
