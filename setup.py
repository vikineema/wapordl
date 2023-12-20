from setuptools import setup, find_packages

setup(
    name = 'wapordl',
    version = '0.3',
    packages = find_packages(include = ['wapordl', 'wapordl.*']),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires = [
        "requests",
        "pandas",
        "numpy",
        "gdal"
    ],
)