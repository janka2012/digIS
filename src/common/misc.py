import errno
import os
import pickle
import re

from os import walk
from definitions import ROOT_DIR


def get_family_names(mypath):
    records = walk(mypath)
    (dirpath, dirnames, filenames) = next(records)
    return dirnames


def get_family_filenames(fam_path, include_subdir=True):
    filename_list = []
    for (dirpath, dirnames, filenames) in walk(fam_path):
        for filename in filenames:
            filename_list.append(os.path.join(dirpath, filename))
        if(not include_subdir):
            break
    return filename_list


def get_filenames_by_extension(dirpath, extension):
    outfiles = []
    list_of_extension = []
    list_of_extension.extend(extension)
    for (dirpath, dirnames, filenames) in os.walk(dirpath):
        for filename in filenames:
            if os.path.splitext(filename)[-1] in list_of_extension:
                outfiles.append(os.path.join(dirpath, filename))
        break

    return outfiles


def get_filenames_by_substring(dirpath, substring):
    outfiles = []
    for (dirpath, dirnames, filenames) in os.walk(dirpath):
        for filename in filenames:
            if substring in filename:
                outfiles.append(os.path.join(dirpath, filename))
        break

    return outfiles


def save_object(obj_name, obj):
    dump_file = os.path.join(ROOT_DIR, "data", "objs", obj_name)
    with open(dump_file, 'wb') as out_handle:
        pickle.dump(obj, out_handle, pickle.HIGHEST_PROTOCOL)


def load_object(obj_name):
    dump_file = os.path.join(ROOT_DIR, "data", "objs", obj_name)
    if os.path.exists(dump_file):
        with open(dump_file, 'rb') as in_handle:
            return pickle.load(in_handle)
    else:
        return None


def delete_file(file):
    try:
        os.remove(file)
    except OSError:
        print("Error while deleting file ", file)


def clean_files(l):
    for i in l:
        if os.path.isfile(i):
            os.remove(i)


def clean_directory(dir_path):
    file_list = os.listdir(dir_path)
    for filename in file_list:
        os.remove(os.path.join(dir_path, filename))


def prepare_directory(dir_path):
    if os.path.exists(dir_path):
        clean_directory(dir_path)
    else:
        os.makedirs(dir_path)


def init_output_dir(output_dir):
    subdirs = ["hmmer", "logs", "pep", "results"]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for subdir in subdirs:
        if not os.path.exists(os.path.join(output_dir, subdir)):
            os.makedirs(os.path.join(output_dir, subdir))


def check_if_file_exists(filename):
    if not os.path.isfile(filename):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
    else:
        return filename


def check_evalue(evalue):

    if evalue is None:
        raise TypeError("Wrong evalue argument type. Current type: {}".format(type(evalue)))

    try:
        evalue_val = float(evalue)

        if evalue_val < 0:
            raise ValueError("Evalue has to be a non-negative number. "
                             "Current value: {}, type: {}".format(evalue, type(evalue)))
    except ValueError:
        raise ValueError("Evalue has to be a non-negative number. "
                         "Current value: {}, type: {}".format(evalue, type(evalue)))


def change_path_to_linux(line):
    matchObj = re.match(r'(.):.*', line, re.M | re.I)
    if matchObj:
        line = line.replace(matchObj.group(1) + ":", '/mnt/' + matchObj.group(1).lower())
    line = line.replace('\\', '/')
    return line
