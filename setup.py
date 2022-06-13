try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit("Please install setuptools.")

import os
import urllib

try:
    from urllib.request import urlretrieve

    classifiers = [
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
except ImportError:
    from urllib import urlretrieve

    classifiers = [
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]

VERSION = "1.0.0"
AUTHOR = "Tyson Dawson, Xinyang Zhang, Ali Rahnavard"
AUTHOR_EMAIL = "gholamali.rahnavard@gmail.com, kathyzhang415@gmail.com"
MAINTAINER = "Tyson Dawson, Xinyang Zhang, Ali Rahnavard"
MAINTAINER_EMAIL = "tysondawson@gmail.com, kathyzhang415@gmail.com, gholamali.rahnavard@gmail.com"

# try to download the bitbucket counter file to count downloads
# this has been added since PyPI has turned off the download stats
# this will be removed when PyPI Warehouse is production as it
# will have download stats
COUNTER_URL = "https://github.com/omicsEye/seqSight/blob/master/README.md"
counter_file = "README.md"
if not os.path.isfile(counter_file):
    print("Downloading counter file to track seqSight downloads" +
          " since the global PyPI download stats are currently turned off.")
    try:
        pass  # file, headers = urlretrieve(COUNTER_URL,counter_file)
    except EnvironmentError:
        print("Unable to download counter")
with open('requirements.txt') as f:
    required = f.read().splitlines()
setup(
    name="seqSight",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    version=VERSION,
    license="MIT",
    description="seqSight: jointly profile microbial strains, genes, and biosynthetic gene clusters from metagenomics data",
    long_description="seqSight: jointly profile microbial strains, genes, and biosynthetic gene clusters from metagenomics data.",
    url="http://github.com/omicsEye/seqSight",
    keywords=['Microbiome', 'metagenomics', 'gene', "biosynthetic gene clusters", "Microbial profiling"],
    platforms=['Linux', 'MacOS', "Windows"],
    classifiers=classifiers,
    # long_description=open('readme.md').read(),
    install_requires=required,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'seqSight = seqSight.seqSight:main'
        ]},
    #test_suite='seqSight.tests.seqSight_test',
    zip_safe=False
)
