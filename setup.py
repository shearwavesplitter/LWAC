import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="LWAC",
    version="0.1.1",
    author="Stefan Mroczek",
    description="Light Weight Auto-Correlator",
    url="https://github.com/shearwavesplitter/LWAC",
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
)
