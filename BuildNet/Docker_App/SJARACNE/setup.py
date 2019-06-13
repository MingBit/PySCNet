from setuptools import setup, find_packages

version = {}
with open("./version.py") as fp:
    exec(fp.read(), version)

setup(
    name="SJARACNe",
    version=version["__version__"],
    description="Gene network reverse engineering from big data",
    url="https://github.com/jyyulab/SJARACNe",
    author="Liang Ding, Alireza Khatamian, Chenxi Qian",
    author_email="liang.ding@stjude.org, chenxi.qian@stjude.org",
    license="See LICENSE.md",
    install_requires=[
        "pandas >= 0.22.0",
        "numpy >= 1.14.2",
        "python-igraph == 0.7.1.post6",
        "scipy == 1.0.1",
        "XlsxWriter == 1.0.2",
    ],
    python_requires="==3.6.1",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    test_suite="tests",
    entry_points={"console_scripts": ["sjaracne=SJARACNe.sjaracne:main"]},
)
