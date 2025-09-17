from os.path import join, abspath, dirname
from setuptools import find_packages, setup


base_dir = abspath(dirname(__file__))
requirements_txt = join(base_dir, 'requirements.txt')
requirements = [l.strip() for l in open(requirements_txt) if l and not l.startswith('#')]

setup(name='eva_release_automation',
      version='0.1.0',
      packages=find_packages(),
      package_data={'run_release_in_embassy': ['nextflow/*'], 'release_automation': ['nextflow/*']},
      install_requires=requirements,
      license='Apache',
      description='EBI EVA - RS Release processing tools',
      url='https://github.com/EBIVariation/eva-release-automation',
      python_requires='>=3.8',
      author='The EVA team',
      author_email='eva-dev@ebi.ac.uk',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Developers',
          'Topic :: Software Development :: Build Tools',
          'License :: OSI Approved :: Apache Software License',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9'
      ]
      )
