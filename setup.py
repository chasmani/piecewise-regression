import pathlib
import setuptools

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setuptools.setup(
    name="piecewise-regression",
    version="0.1.1",
    description="piecewise (segmented) regression in python",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/chasmani/piecewise-regression",
    author="Charlie Pilgrim",
    author_email="pilgrimcharlie2@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=setuptools.find_packages(),
    install_requires=["numpy", "matplotlib", "scipy", "statsmodels"],
)