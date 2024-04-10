from setuptools import setup, find_packages

setup(
    name = 'wapordl',
    version = '0.11.0',
    packages = find_packages(include = ['wapordl', 'wapordl.*']),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires = [
        "requests",
        "pandas>=2.1.0,<3",
        "numpy>=1.15,<2",
        "gdal>=3.4.0,<4",
        "shapely>=2.0.0",
        "tqdm",
    ],
)