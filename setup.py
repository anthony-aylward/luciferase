import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='luciferase',
    version='1.8',
    author='Anthony Aylward',
    author_email='aaylward@eng.ucsd.edu',
    description='Helper functions for luciferase data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/anthony-aylward/luciferase.git',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    install_requires=['matplotlib', 'pandas', 'scipy', 'seaborn'],
    entry_points={
        'console_scripts': ['reporter-barplot=luciferase.luciferase:main',]
    }
)
