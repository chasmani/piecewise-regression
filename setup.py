import setuptools

# This call to setup() does all the work
setuptools.setup(
    name="piecewise-regression",
    version="1.1.3",
    description="piecewise (segmented) regression in python",
    long_description=   "piecewise-regression provides tools for fitting " 
                        "continuous straight line models to data with "
                        "breakpoint(s) where the gradient changes. "
                        ""
                        "For docs and more information, "
                        "visit the Github repo at "
                        "https://github.com/chasmani/piecewise-regression.",
    long_description_content_type="text/markdown",
    url="https://github.com/chasmani/piecewise-regression",
    author="Charlie Pilgrim",
    author_email="pilgrimcharlie2@gmail.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    packages=setuptools.find_packages(),
    install_requires=["numpy", "matplotlib", "scipy", "statsmodels"],
)

