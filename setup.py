from __future__ import print_function
from setuptools import setup, find_packages
from setuptools.command.install import install
import sysconfig, os, shutil, sys
import subprocess as sp

# based on: https://stackoverflow.com/questions/431684/how-do-i-change-directory-cd-in-python
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, new_path):
        self.new_path = os.path.expanduser(new_path)

    def __enter__(self):
        self.saved_path = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.saved_path)

# based on:
# * https://stackoverflow.com/questions/33168482/compiling-installing-c-executable-using-pythons-setuptools-setup-py
# * https://blog.niteo.co/setuptools-run-custom-code-in-setup-py/
def setup_compiler():
    for var in ['CC', 'CXX', 'CPPFLAGS', 'CFLAGS', 'LDFLAGS']:
        if var in os.environ:
            print('# rfmix: (env) {}={}'.format(var, os.environ[var]))
        elif var in sysconfig.get_config_vars():
            value = sysconfig.get_config_var(var)
            print("# rfmix: (sysconfig) {}={}".format(var, value))
            os.environ[var] = value

def get_virtualenv_path():
    """Used to work out path to install compiled binaries to."""
    if hasattr(sys, 'real_prefix'):
        return sys.prefix

    if hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        return sys.prefix

    if 'conda' in sys.prefix:
       return sys.prefix

    return None

def compile_pop_phased():
    src_path = './PopPhased/'
    cxx = os.environ['CXX'] # c++ compiler
    cxx_options = ['-Wall', '-O3', '-ftree-vectorize', '-fopenmp']
    srcs = [
        'main.cpp',
        'getdata.cpp',
        'randomforest.cpp',
        'crfviterbi.cpp',
        'windowtosnp.cpp'
    ]
    output = 'RFMix_PopPhased'
    cmd = "{cxx} {options} {srcs} -o {out}".format(
        cxx=cxx,
        options=' '.join(cxx_options),
        srcs=' '.join(srcs),
        out=output
    )
    dst_path = os.path.join(sys.prefix, 'bin')
    with cd(src_path):
        print("Inside {}".format(os.getcwd()))
        print('# rfmix: (pop-phased cmd) {}'.format(cmd))
        sp.check_call(cmd, shell=True)
        print('# rfmix: (pop-phased) {} -> {}'.format(output, dst_path))
        shutil.copy2(output, dst_path)

def compile_trio_phased():
    src_path = './TrioPhased/'
    cxx = os.environ['CXX'] # c++ compiler
    cxx_options = ['-Wall', '-O3', '-ftree-vectorize', '-fopenmp']
    srcs = [
        'main.cpp',
        'getdata.cpp',
        'randomforest.cpp',
        'crfviterbi.cpp',
        'windowtosnp.cpp'
    ]
    output = 'RFMix_TrioPhased'
    cmd = "{cxx} {options} {srcs} -o {out}".format(
        cxx=cxx,
        options=' '.join(cxx_options),
        srcs=' '.join(srcs),
        out=output
    )
    dst_path = os.path.join(sys.prefix, 'bin')
    with cd(src_path):
        print("Inside {}".format(os.getcwd()))
        print('# rfmix: (trio-phased cmd) {}'.format(cmd))
        sp.check_call(cmd, shell=True)
        print('# rfmix: (tro-phased) {} -> {}'.format(output, dst_path))
        shutil.copy2(output, dst_path)

def compile_and_install_software():
    """Use the subprocess module to compile/install the C software."""
    setup_compiler()
    compile_pop_phased()
    compile_trio_phased()

class CustomInstall(install):
    def run(self):
        compile_and_install_software()
        install.run(self)

setup(
    name="RFMix",
    version="1.5.4a",
    packages=find_packages(exclude=('TestData', 'docs')),
    scripts=['RunRFMix.py'],
    author='modified: Hall Lab, Washington University in St. Louis, original software: Bustamante Lab, Stanford University',
    author_email='idas@wustl.edu',
    description='A discriminative method for local ancestry inference',
    url='https://sites.google.com/site/rfmixlocalancestryinference/',
    python_requires='>=2.6, <3',
    include_package_data=True,
    cmdclass={'install': CustomInstall}
)
