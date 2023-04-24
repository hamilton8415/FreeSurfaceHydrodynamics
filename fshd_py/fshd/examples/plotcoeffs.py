#!/usr/bin/pyton3

from fshd import FS_HydroDynamics


def main():
    import importlib.resources as ilr
    import subprocess
    import time

    rho = 1025.
    g = 9.81
    buoy_mass = 1400.  # kg
    BuoyA5 = FS_HydroDynamics(1.0, g, rho)

    dt = 0.005

    # Set area and 2nd moments of area for waterplane
    BuoyA5.SetWaterplane(5.47, 1.37, 1.37)
    # Set COB relative to waterplane coordinate system.
    BuoyA5.SetCOB(0., 0., -0.22)
    #Set COG relative to waterplane coordinate system.
    BuoyA5.SetCOG(0., 0., -0.24)
    BuoyA5.SetVolume(buoy_mass / rho)
    BuoyA5.SetMass(buoy_mass)

    HydrodynamicsBaseFilename = \
        ilr.files('fshd.examples').joinpath('example_hydrodynamic_coeffs',
                                            'BuoyA5')
    HydrodynamicsBaseFilename = str(HydrodynamicsBaseFilename)
    BuoyA5.ReadWAMITData_FD(HydrodynamicsBaseFilename)
    BuoyA5.ReadWAMITData_TD(HydrodynamicsBaseFilename)
    BuoyA5.SetTimestepSize(dt)

    BuoyA5.Plot_FD_Coeffs()
    BuoyA5.Plot_TD_Coeffs()

    try:
        print("Enter Ctrl-C to quit. (Enter 'pkill gnuplot_qt' to clear plots if necessary)")
        while True:
            time.sleep(0.1)
    except KeyboardInterrupt:
        subprocess.run(['pkill', 'gnuplot_qt'])


if __name__ == '__main__':
    main()
