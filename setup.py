'''
Created on May 15, 2014

@author: t.dengg
'''

from distutils.core import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='pylastic',
      version='0.1',
      description='Calculate ab-initio elastic constants.',
      url='http://github.com/tdengg/pylastic',
      author='Thomas Dengg',
      author_email='thomas.dengg@mcl.at',
      license='LGPL',
      packages=['pylastic'],
      install_requires=['numpy','lxml','matplotlib','json'],
      zip_safe=False)
