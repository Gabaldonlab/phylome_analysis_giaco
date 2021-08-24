from setuptools import setup, find_packages


with open("README.md") as f:
    readme = f.read()

with open("LICENSE.txt") as f:
    license = f.read()

setup(
    name="analyse_phylome",
    version="0.0.1",
    description="Analyse phylome data",
    long_description=readme,
    url="https://github.com/giacomomutti/analyse_phylome",
    author="Giacomo Mutti",
    author_email="giaco98@hotmail.it",
    license=license,
    packages=find_packages(),  # exclude=('tests', 'docs')),
    zip_safe=False,
    include_package_data=True,
)
