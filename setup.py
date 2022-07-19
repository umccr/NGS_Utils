#!/usr/bin/env python
import os
from os.path import join
from setuptools import setup

pkg = 'ngs_utils'

try:
    import versionpy
except ImportError:
    res = input('Installation requires versionpy. Install it now? [Y/n]')
    if res.lower().startswith('n'):
        raise
    os.system('pip install versionpy')
    import versionpy


package_data = {
    pkg: versionpy.find_package_files('', pkg, skip_exts=['.sass', '.coffee'])
}


setup(
    name=pkg,
    version='2.9.2',
    author='Vlad Savelyev',
    author_email='vladislav.sav@gmail.com',
    description='Python utilities for bioinformatics tools and pipelines',
    long_description=open('README.md').read(),
    url='https://github.com/umccr/NGS_Utils',
    license='GPLv3',
    packages=[pkg],
    package_data=package_data,
    include_package_data=True,
    zip_safe=False,
    install_requires=['click', 'versionpy'],
    scripts=[path for path in
         [join('scripts', fn) for fn in [
             'standardize_bed.py',
             'split_bed.py',
             'qstat.py',
             'sort_bed.py',
             'bed_file_from_gene_list.py',
             'html_to_base64.py',
             'group_and_merge_by_gene.py',
             'hg19_addchr.py',
             'generate_bed.py',
         ]]],
    keywords='bioinformatics',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
