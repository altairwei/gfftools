from setuptools import setup, find_packages

setup(
    name="GFFTool",
    version="0.0.1",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "HTSeq"
    ],
    entry_points='''
        [console_scripts]
        gfftool=gfftool.main:cli
    ''',
)