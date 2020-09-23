Installation
===================

Requirements
--------------
.. note::
    As pyscnet integrates `docker`_ for gene regulatory construction, it is necessary to install docker before the installation. Check here for `docker installation`_.


PyPI
--------------
install pyscnet via PyPI

.. code-block::

    pip install pyscnet

.. warning::

    It might occurs package version conflicts.

To address this issue, you can create a new/clean conda environment

.. code-block::

    #create new environment with python3.6
    conda create --name pyscnet_env python=3.6

    #activate pyscnet_env
    conda activate pyscnet_env

    #install pyscnet
    pip install pyscnet

    #go to python
    import pyscnet


Github
--------------
install develop version pyscnet via github

.. code-block::

    git clone https://github.com/MingBit/PySCNet
    mkdir dist | python setup.py sdist
    pip install dist/pyscnet-0.0.3.tar.gz


.. _docker: https://www.docker.com/
.. _docker installation: https://docs.docker.com/get-docker/