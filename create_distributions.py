"""Script to automate the creation of Windows and Linux source distributions.

The TOPKAPI_example directory is also copied and the .svn directories stripped
to make a clean distribution. The manual is included in MSWord format for now
because this is how it's stored in SVN.

This script currently relies on Linux tools and will only work on a Linux
system.

"""
import os
import shutil

# Linux shell command to strip .svn directories
command = 'find . -name .svn -type d -print0 | xargs -0 rm -rf'

def make_distro(dist_path, ex_path, add_files):
    path = os.path.join(dist_path, ex_path)

    if os.path.isdir(dist_path):
        for root, dirs, files in os.walk(dist_path, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
    else:
        os.mkdir(dist_path)

    shutil.copytree(ex_path, path)
    curr_dir = os.getcwd()
    os.chdir(path)
    os.system(command)
    os.chdir(curr_dir)

    for fname in add_files:
        shutil.copy(fname, dist_path)

if __name__ == "__main__":
    # build source distributions
    os.system('python setup.py sdist --formats=gztar,zip --force-manifest')

    # make Linux distribution
    dist_path = 'TOPKAPI_linux'
    ex_path = 'TOPKAPI_Example'
    linux_files = ['dist/TOPKAPI-0.1.1.tar.gz', 'TOPKAPI_Manual.doc']

    make_distro(dist_path, ex_path, linux_files)

    # make Windows distribution
    dist_path = 'TOPKAPI_windows'
    ex_path = 'TOPKAPI_Example'
    windows_files = ['dist/TOPKAPI-0.1.1.zip', 'TOPKAPI_Manual.doc']

    make_distro(dist_path, ex_path, windows_files)
