import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='luciferase',
    version='1.5.5',
    author='Anthony Aylward, Joshua Chiou, Mei-Lin Okino',
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
    install_requires=[
        'estimateratio', 'matplotlib', 'pandas', 'scipy', 'seaborn'
    ],
    entry_points={
        'console_scripts': [
            'luciferase-barplot=luciferase.luciferase:main',
            'luciferase-ratioplot=luciferase.ratioplot:main'
        ]
    }
)
