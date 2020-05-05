"""place for file handling functions"""
# -*- coding: utf-8 -*-


import os
import shutil
from nephele2.pipelines.pipebase import PipeBase


def cat_file(cat_from, cat_to):
    """cat_file

    appends one file to another.
    :param cat_from:
    :param cat_to:
    """
    with open(cat_to, 'a') as _a_file:
        with open(cat_from, 'r') as _r_file:
            shutil.copyfileobj(_r_file, _a_file)


def rm_sym_links_from(dname):
    """
    Removes all symbolic links from the specified directory.

    Args:
        dname (str): the path of the directory

    """
    for fname in os.listdir(dname):
        filepath = os.path.join(dname, fname)
        if os.path.islink(filepath):
            os.unlink(filepath)
    return os.listdir(dname)


def remove_intermediate_dirs(dname, dirs_to_keep):
    """remove_intermediate_dirs
    Removes dirs_to_keep from dname
    :param dname:
    :param dirs_to_keep:
    """
    output_dirs = os.listdir(dname)
    for output_dir in output_dirs:
        if output_dir in dirs_to_keep:
            continue
        dir_path = os.path.join(dname, output_dir)
        if os.path.isdir(dir_path):
            shutil.rmtree(dir_path, ignore_errors=True)


def remove_intermediate_files(dname, files_to_keep):
    """Remove files not in `FILES_TO_KEEP` or `DIRS_TO_KEEP`.

    Args:
      dname (str): path to output directory
      map_fp (str): map filename

    """
    output_fnames = rm_sym_links_from(dname)
    for output_fname in output_fnames:
        if output_fname in files_to_keep:
            continue
        output_fp = os.path.join(dname, output_fname)
        if os.path.isfile(output_fp):
            os.remove(output_fp)


def get_current_files(dname):
    """
    Gets the current files saved by mothur from the current_files.summary file
    generated by mothur's get.current() command.

    Args:
        dname (str): the path of the outputs directory

    Returns:
        dict(str,str): a dictionary of the current file set::

            {
                'file_type1': file_path1,
                'file_type2': file_path2
            }
    """
    current_summ = dname + 'current_files.summary'
    PipeBase.ensure_file_exists(current_summ)
    current_set = {}
    for line in open(current_summ, "r"):
        ftype, fpath = line.split("=")
        current_set[ftype.strip()] = fpath.strip()
    return current_set