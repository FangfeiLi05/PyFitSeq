from setuptools import setup

with open('README.md','r') as fh:
    long_description = fh.read()

setup(
    name='PyFitSeq',
    version='1.2.0-rc1',
    description='python version of FitSeq (the latest)',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=['PyFitSeq'],
    scripts=['PyFitSeq/evo_simulator.py','PyFitSeq/pyfitseq.py'],
    install_requires=['numpy>=1.17.3','pandas>=0.25.3','scipy>=1.3.1'],
    license='MIT',
    long_description=long_description,
    long_description_content_type='text/markdown'
    )
