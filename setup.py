from setuptools import setup

exec(open('main/version.py').read())

setup(
    name='CARTtools',
    version=__version__,
    description='CARTtools is ...',
    url='x@x',
    author='RahmanTeam',
    author_email='rahmanlab@icr.ac.uk',
    licence='MIT',
    packages=['main', 'cart2enst_'],
    scripts=['bin/CART2ENST.py'],
    zip_safe=False
)
