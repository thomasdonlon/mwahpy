import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mwahpy",
    version="1.3.4",
    author="Tom Donlon",
    author_email="donlot@rpi.edu",
    description="A python package for easily parsing and processing data from MilkyWay@home",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/thomasdonlon/mwahpy",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.2.3',
    install_requires=[
        "numpy",
        "galpy",
        "matplotlib",
        "astropy",
        "scipy"
        ],
)
