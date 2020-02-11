import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biomassmw", # Replace with your own username
    version="0.0.1",
    author="Joshua Chan",
    author_email="joshua.chan@example.com",
    description="Package for determining the chemical formulae and molecular weights of macromolecules in genome-scale metabolic models",
    long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/chan_csu/biomassMWpy",
    packages=["biomassmw", "biomassmw.objects"],    #setuptools.find_packages()],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
