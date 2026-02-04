from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="recmpox",
    version="0.0.1",
    author="RecMpox contributors",
    description="RecMpox: Classify consensus mpox genomes at diagnostic SNPs (recombinant calling).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[],  # matplotlib optional (for --plot)
    entry_points={
        "console_scripts": [
            "recmpox=recmpox.recmpox:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
