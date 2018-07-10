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
    packages=[
        'main',
        'select_nm_',
        'map_nm_',
        'select_enst_',
        'format_cart_',
        'pipeline'
    ],
    scripts=[
        'bin/select_NM.py',
        'bin/map_NM.py',
        'bin/select_ENST.py',
        'bin/format_CART.py',
        'bin/CART_pipeline.py'
    ],
    zip_safe=False
)
