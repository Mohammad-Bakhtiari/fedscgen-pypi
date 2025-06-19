from setuptools import setup, find_packages

setup(
    name="fedscgen",
    version="0.0.4",
    packages=find_packages(),
    install_requires=[
        "anndata==0.9.1",
        "matplotlib==3.7.2",
        "numpy==1.24.4",
        "pandas==1.5.3",
        "scanpy==1.9.3",
        "scipy==1.11.1",
        "seaborn==0.12.2",
        "scikit-learn==1.3.0",
        "h5py==3.9.0",
        "pyyaml==6.0",
        "requests==2.31.0",
        "torch==2.7.0",
        "torchvision==0.22.0",
        "pytorch-lightning==2.0.5",
        "scgen==2.1.0",
        "pyro-ppl==1.8.5",
        "python-louvain==0.16",
        "rpy2==3.5.13",
    ],
    entry_points={
        "console_scripts": [
            "fedscgen-centralized=fedscgen.cli.run_centralized:main",
            "fedscgen-federated=fedscgen.cli.run_federated:main",
        ]
    },
    author="Mohammad Bakhtiari",
    author_email="mohammad.bakhtiari@uni-hamburg.de",
    description="A package for single-cell genomics analysis with scGen",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/Mohammad-Bakhtiari/fedscgen-pypi",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.10,<3.11",
)