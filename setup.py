from setuptools import setup, find_packages

setup(
    name = 'wapordl',
    version = '0.9',
    packages = find_packages(include = ['wapordl', 'wapordl.*']),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires = [
        "requests",
        "pandas>=2.1.0",
        "numpy",
        "gdal>=3.4.0",
        "shapely>=2.0.0",
        "tqdm",
    ],
)