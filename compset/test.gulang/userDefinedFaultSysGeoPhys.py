#! /usr/env/bin python3
def defineSysPhys(ftSystem, ftNames, xCoorDict, yCoorDict):
    ftPhys = {}
    for ftNameKey in ftNames:
        tmp = []
        for i, xcoor in enumerate(xCoorDict[ftNameKey]):

            if ftSystem=="none":
                ftType = 1 # 1: left-strike; -1: right-strike; 2: thrust; -2: normal.
                ftDip = 90 # 90: strike-slip; positive: tilting to y+; negative; tilting to y-.
                ftLoadMaxShear = 1.427e-14
                ftLoadAngle = -999
                ftLoadWt = 1.
                ftVis = 6e21
                tmp += [[ftType, ftDip, ftLoadMaxShear, ftLoadAngle, ftLoadWt, ftVis]]
            elif ftSystem=="gulang":
                ftType = 1 # left-lateral strike-slip
                ftDip = 90
                ftLoadMaxShear = 1.427e-14
                ftLoadAngle = -999
                ftLoadWt = 1.
                ftVis = 6e21

                tmp += [[ftType, ftDip, ftLoadMaxShear, ftLoadAngle, ftLoadWt, ftVis]]

        ftPhys[ftNameKey] = tmp
    return ftPhys

def plotSystemPhys(ftPhys, ftNames, xCoorDict):
    fig, ax = plt.subplots(3, 2, figsize=(15, 10), dpi=600)

    for key in ftNames:
        ax[0,0].plot(xCoorDict[key], [row[0] for row in ftPhys[key]])
        ax[0,1].plot(xCoorDict[key], [row[1] for row in ftPhys[key]])
        ax[1,0].plot(xCoorDict[key], [row[2] for row in ftPhys[key]])
        ax[1,1].plot(xCoorDict[key], [row[3] for row in ftPhys[key]])
        ax[2,0].plot(xCoorDict[key], [row[4] for row in ftPhys[key]])
        ax[2,1].plot(xCoorDict[key], [row[5] for row in ftPhys[key]])
    #plt.savefig('mesh.png', dpi=600)
    plt.show()
