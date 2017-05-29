import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='scRNApipe',
    version='0.1.0',
    description='Package for analysing scRNA-seq in Transcript Tag Counting data.',
    long_description=read('README.md'),
    author='Stavros Giannoukakos',
    author_email='s.p.giannoukakos@hotmail.com',
    packages=['scRNApipe'],
    url=['https://github.com/MarinusVL/scRNApipe'],
    keywords=['single cell RNA analysis'],

    install_requires=['pysam>=0.8.3', 'numpy', 'multiqc', 'STAR', 'umis', 'umi_tools', ,'python>=2.5,<3','natsort'],

    dependency_links=['https://sourceforge.net/projects/subread/files/subread-1.5.2/subread-1.5.2-source.tar.gz/download',
                      'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5_source.zip'
                      ],
    package_data = {
        '': ['configuration_file.txt']
        },
    entry_points={
          'console_scripts': ['scRNApipe = scRNApipe.scRNApipe:main']
      },
)
