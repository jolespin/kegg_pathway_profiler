from setuptools import setup, find_packages

from os import path

script_directory = path.abspath(path.dirname(__file__))

package_name = "kegg_pathway_profiler"
version = None
with open(path.join(script_directory, package_name, '__init__.py')) as f:
    for line in f.readlines():
        line = line.strip()
        if line.startswith("__version__"):
            version = line.split("=")[-1].strip().strip('"')
assert version is not None, f"Check version in {package_name}/__init__.py"

requirements = list()
with open(path.join(script_directory, 'requirements.txt')) as f:
    for line in f.readlines():
        line = line.strip()
        if line:
            if not line.startswith("#"):
                requirements.append(line)

from setuptools import setup, find_packages

setup(
    name='kegg_pathway_profiler',
    version=version,
    description='KEGG Pathway Profiler',
    long_description=open(path.join(script_directory, 'README.md')).read(),
    long_description_content_type='text/markdown',
    author='Josh L. Espinoza',
    url='https://github.com/jolespin/kegg_pathway_profiler',
    license='GPLv3',
    packages=find_packages(),  # Automatically find packages in the project
    include_package_data=True,  # Include package data as defined in MANIFEST.in
    # package_data={
    #     'kegg_pathway_profiler': ['data/database.pkl.gz'],  # Specify the data file to include
    # },
    install_requires=[
    ],
    scripts=[
        "bin/profile-pathway-coverage.py",
        "bin/build-pathway-database.py",
        "bin/download-kegg-pathways.sh",
    ],
    requirements=requirements,
    python_requires='>=3.6',
)