import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = ["cobra>=0.15.1"]

setuptools.setup(
    name="biomassmw",
    version="1.0.0",
    author="Joshua Chan",
    author_email="joshua.chan@colostate.edu",
    description="Package for determining the chemical formulae and molecular weights of macromolecules in genome-scale metabolic models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chan_csu/biomassMWpy",
    packages=["biomassmw", "biomassmw.objects"],
    install_requires=requirements,
    classifiers=[
        "Operating System :: OS Independent",
    ]

)
