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

VERSION = "1.0"
AUTHOR = "Ali Rahnavard, Xinyang Zhang"
AUTHOR_EMAIL = "gholamali.rahnavard@gmail.com, kathyzhang415@gmail.com"
MAINTAINER = "Ali Rahnavard, Xinyang Zhang"
MAINTAINER_EMAIL = "gholamali.rahnavard@gmail.com, kathyzhang415@gmail.com"

# try to download the bitbucket counter file to count downloads
# this has been added since PyPI has turned off the download stats
# this will be removed when PyPI Warehouse is production as it
# will have download stats
COUNTER_URL = "https://github.com/omicsEye/deepStrain/blob/master/README.md"
counter_file = "README.md"
if not os.path.isfile(counter_file):
    print("Downloading counter file to track seqSight downloads" +
          " since the global PyPI download stats are currently turned off.")
    try:
        pass  # file, headers = urlretrieve(COUNTER_URL,counter_file)
    except EnvironmentError:
        print("Unable to download counter")

setup(
    name="seqSight",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    version=VERSION,
    license="MIT",
    description="seqSight: multi-resolution clustering",
    long_description="seqSight provides an elegant clustering approach " + \
                     "to find clusters in data sets with different density and resolution.",
    url="http://github.com/omicsEye/deepStrain",
    keywords=['clustering', 'multi-resolution', 'hierarchically'],
    platforms=['Linux', 'MacOS', "Windows"],
    classifiers=classifiers,
    # long_description=open('readme.md').read(),
    install_requires=[
        "latex >= 0.0.1",
        "Cython >= 0.29.2",
        "Numpy >= 1.9.2",
        "Scipy >= 0.17.0",
        "Matplotlib >= 1.5.1",
        "Scikit-learn >= 0.14.1",
        "pandas >= 0.18.1",
        "Community == 1.0.0b1",
        "networkx >= 0.2.5"
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'seqSight = seqSight.seqSight:main',
            'deepStrainviz = seqSight.viz:main',
            'deepStrain_test = seqSight.tests.deepStrain_test:main'
        ]},
    test_suite='seqSight.tests.deepStrain_test',
    zip_safe=False
)
