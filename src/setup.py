import setuptools

setuptools.setup(
    name='vsat',
    version='1.0',
    scripts=['./scripts/vsat'],
    author='Kijin Kim',
    description='Virus Sequence Analysis Tools',
    packages=['vsat'],
    install_requires=[
        'setuptools',
        'requests',
        'biopython',
        'wget',
    ],
    include_package_data=True
)
