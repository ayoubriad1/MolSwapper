from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="MolSwapper",
    version="1.0.0",
    author="Ayoub",
    author_email="t61199988@gmail.com",
    description="A modular toolkit for molecular fragment swapping and chimeric structure generation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ayoubriad1/MolSwapper",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.9",
    install_requires=[
        "rdkit>=2023.09.1",
        "Pillow>=10.0.0",
    ],
    entry_points={
        "console_scripts": [
            "molswapper=molswapper.main:main",
        ],
    },
)