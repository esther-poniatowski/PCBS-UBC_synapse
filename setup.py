#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Makes the projet installable with : `pip install -e`."""

import setuptools

with open("README.md", "r") as long_descr :
    long_description = long_descr.read()

setuptools.setup(
	name="UBC_synapse",
	version="0.0.1",
	author="esther-poniatowski",
	author_email="esther.poniatowski@ens.fr",
	description="",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/esther-poniatowski/UBC_synapse",
	packages=setuptools.find_packages("UBC_synapse"),
	package_dir={"": "UBC_synapse"},
	package_data={"": "UBC_data/*"},
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	python_requires='>=3.6',
	install_requires=['matplotlib',
					'numpy',
					'pandas',
					'os',
					'ctypes',
					'subprocess'
					'pickle',
					'warnings'])
)