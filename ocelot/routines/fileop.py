import shutil
import os
import glob
import re
import pickle


def write_obj(obj, fname):
    """
    write an object into binary with pickle

    :param obj:
    :param fname:
    :return:
    """
    with open(fname, 'wb') as f:
        pickle.dump(obj, f)


def read_obj(fname):
    """
    read an object from pickle

    :param fname:
    :return:
    """
    with open(fname, 'rb') as f:
        obj = pickle.load(f)
    return obj


def check_path(pathlist):
    """
    check if the paths in the pathlist are all valid

    :param pathlist:
    :return:
    """
    normal = 1
    for i in pathlist:
        if not os.path.isdir(i):
            normal = 0
    return normal


def movefile(what, where):
    """
    shutil operation to move

    :param what:
    :param where:
    :return:
    """
    try:
        shutil.move(what, where)
    except IOError:
        os.chmod(where, 777)
        shutil.move(what, where)

def removefile(what):
    try:
        os.remove(what)
    except OSError:
        pass


def copyfile(what, where):
    """
    shutil operation to copy

    :param what:
    :param where:
    :return:
    """
    try:
        shutil.copy(what, where)
    except IOError:
        os.chmod(where, 777)
        shutil.copy(what, where)


def createdir(directory):
    """
    mkdir

    :param directory:
    :return:
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def lsexternsion(path, extension='cif'):
    """
    get a list of file names with certain extension in the path

    :param path:
    :param extension:
    :return: [path/filename1.extension, ...]
    """
    fns = glob.glob(path + '/*.' + extension)
    return fns


def nonblank_lines(f):
    """
    get nonblank lines in a file

    :param f: string, file name
    :return:
    """
    for l in f:
        line = l.rstrip()
        if line:
            yield line


def check_outcar(filename='OUTCAR'):
    """
    check if OUTCAR normally term

    :param filename:
    :return:
    """
    normal = 0
    try:
        outcar = open(filename, 'r')
    except IOError:
        return normal
    for line in outcar:
        if re.search('Total CPU time used ', line):
            normal = 1
    return normal


def check_gausslog(filename):
    """
    check if gaussian log normally term

    :param filename:
    :return:
    """
    normal = 0
    try:
        out = open(filename, 'r')
    except IOError:
        return normal
    for line in out:
        if re.search('Normal termination', line):
            normal = 1
    return normal


def check_allgausslog(howmanylog):
    """
    check if certain amount of gaussian job finished

    :param howmanylog:
    :return:
    """
    normal = 0
    finfile = 0
    fns = glob.glob('*.log')
    if len(fns) == 0:
        return normal
    for filename in fns:
        try:
            out = open(filename, 'r')
        except IOError:
            return normal
        for line in out:
            if re.search('Normal termination', line):
                finfile += 1
    if finfile == howmanylog:
        normal = 1
    return normal

import inspect

def retrieve_name(var):
    """
    Gets the name of var. Does it from the out most frame inner-wards.
    https://stackoverflow.com/questions/18425225/getting-the-name-of-a-variable-as-a-string

    :param var: variable to get name from.
    :return: string
    """
    for fi in reversed(inspect.stack()):
        names = [var_name for var_name, var_val in fi.frame.f_locals.items() if var_val is var]
        if len(names) > 0:
            return names[0]