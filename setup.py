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
        'selectnms_',
        'mapnms_',
        'selectensts_',
        'formatcarts_',
        'pipeline',
        'refseqdb_',
        'refseqcheck_',
        'ensembldb_',
        'compareensts_',
        'summarize_'
    ],
    scripts=[
        'bin/SelectNMs.py',
        'bin/MapNMs.py',
        'bin/SelectENSTs.py',
        'bin/FormatCARTs.py',
        'bin/CARTpipeline.py',
        'bin/RefSeqDB.py',
        'bin/RefSeqCheck.py',
        'bin/EnsemblDB.py',
        'bin/CompareENSTs.py',
        'bin/Summarize.py'
    ],
    zip_safe=False
)
