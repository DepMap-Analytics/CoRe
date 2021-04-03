import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    requirements = f.read().split("\n")

setuptools.setup(
    name="CoRe",
    version="0.1.0",
    author="Emre Karakoc, Clare Pacini, Alessandro Vinceti and Francesco Iorio",
    author_email="francesco.iorio@sanger.ac.uk",
    description="CoRe",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DepMap-Analytics/CoRe",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

