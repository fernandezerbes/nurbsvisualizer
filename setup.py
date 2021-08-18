from setuptools import setup

setup(name='nurbsvisualizer',
      version='0.1',
      description='Visualization of rational and non-rational B-Splines',
      url='https://github.com/FernandezErbes/nurbsvisualizer',
      author='Federico Fernández Erbes',
      author_email='fernandezerbes@gmail.com',
      license='MIT',
      packages=['nurbsvisualizer'],
      python_requires='>=3',
      install_requires=[
        'numpy>=1.18.1',
        'matplotlib>=3.2.2'],
      zip_safe=False)
