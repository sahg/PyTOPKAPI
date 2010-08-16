from distutils.core import setup

setup(name='TOPKAPI',
      version='0.2dev',
      description='SAHG TOPKAPI model implementation',
      author='Theo Vischel & Scott Sinclair',
      author_email='theo.vischel@hmg.inpg.fr; sinclaird@ukzn.ac.za',
      packages=['TOPKAPI',
                'TOPKAPI.parameter_utils',
                'TOPKAPI.results_analysis'],
      package_dir={'TOPKAPI':'lib'}
      )

