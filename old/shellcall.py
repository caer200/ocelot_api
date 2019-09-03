import os
import warnings
import subprocess
import fileop


def aflow(aflow_p, jobtype, poscar='POSCAR'):
    """
    call aflow to std/convert POSCAR
    csplit is used for std
    :param aflow_p:
    :param jobtype:
    :param poscar:
    :return:
    """
    if jobtype == 'frac2kart':
        out = ['cart.poscar']
        os.system(aflow_p + ' --cart <' + poscar + '>' + out)
        return out
    elif jobtype == 'kart2frac':
        out = ['frac.poscar']
        os.system(aflow_p + ' --frac <' + poscar + '>' + out)
        return out
    elif jobtype == 'std':
        os.system(aflow_p + " --kpath < " + poscar + " > aflow_kpath")
        os.system("csplit --silent --prefix=kpath aflow_kpath ''\'/^\/\/\s/1\''' ''\'{2}\''' ")
        os.system("rm kpath00 kpath03 aflow_kpath")
        os.system("sed -i ''\'$ d\''' kpath* ")
        fileop.movefile('kpath01', 'std_unwrap.poscar')
        fileop.movefile('kpath02', 'kpath_unwrap.poscar')
        out = ['std_unwrap.poscar', 'kpath_unwrap.poscar']
        return out
    else:
        warnings.warn('jobtype in function SysCall.aflow is not recognized')


def oszicar_energy(filename='OSZICAR'):
    """
    get energy from OSZICAR, using tail
    :param filename:
    :return:
    """
    line = subprocess.check_output(['tail', '-1', filename])
    entry = str(line.decode('unicode_escape')).split()
    en_str = entry[2]
    try:
        energy = float(en_str)
    except ValueError:
        # print ('not finished at ' + filename)
        energy = 0.0
    return energy


def runvasp(vasp='vasp_std'):
    """
    run vasp, using mpirun, make sure correct modules are loaded
    :param vasp:
    :return:
    """
    if not os.path.isfile('OUTCAR'):
        os.system('mpirun ' + vasp)
    else:
        warnings.warn('there is already an OUTCAR at ' + os.getcwd())


def rungauss(comfile, gauss='g16'):
    """
    run gaussian, technically .gjf can also be used
    :param comfile:
    :param gauss:
    :return:
    """
    os.system(gauss + ' ' + comfile)
    # runJob = subprocess.Popen([gauss, comfile])
    # runJob.wait()
