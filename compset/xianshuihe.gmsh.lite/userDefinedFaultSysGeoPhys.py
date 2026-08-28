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
                ftLoadWt = 450.
                ftVis = 8.4e21
                tmp += [[ftType, ftDip, ftLoadMaxShear, ftLoadAngle, ftLoadWt, ftVis]]
            elif ftSystem=="xianshuihe":
                ftType = 1   # left-lateral strike-slip (Xianshuihe)
                ftDip = 90   # all faults vertical
                # Qiao et al. (2022) Fig. 4b give a SLIP rate of ~5-13 mm/yr,
                # not a strain rate. ftLoadMaxShear sets an asymptotic stress
                # in interstress.f90, so the long-term slip rate is emergent
                # and has to be calibrated. Starting value scales the SAF
                # anchor (1.427e-14 -> ~38 mm/yr peak) to a ~13 mm/yr target:
                #   1.427e-14 * 13/38 = 4.9e-15
                # Then measure with plot_saf_figure3.py and rescale linearly.
                ftLoadMaxShear = 4.9e-15
                ftLoadAngle = -999
                ftLoadWt = 450.
                ftVis = 5.0e21   # starting value, as tuned for gulang
                tmp += [[ftType, ftDip, ftLoadMaxShear, ftLoadAngle, ftLoadWt, ftVis]]
            elif ftSystem=="gulang":
                ftType = 1 # left-lateral strike-slip
                ftDip = 90
                ftLoadMaxShear = 1.427e-14
                ftLoadAngle = -999
                ftLoadWt = 450.
                ftVis = 5.0e21   # 1e21 was too low (no nucleation), 8.4e21 too high
                                 # (tensile normal at large strike). 5e21 keeps
                                 # γ·ant·sin(2φ) below 100 MPa ambient AND gets
                                 # rs∞ above 50 MPa nucleation threshold.

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
