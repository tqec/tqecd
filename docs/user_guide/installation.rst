How to install ``tqecd``
=======================

Requirements before installing
------------------------------

Python version
~~~~~~~~~~~~~~

The ``tqecd`` package only supports Python 3.10 and onward. If you have Python 3.9 or below,
please update your Python installation.

Additional toolchains
~~~~~~~~~~~~~~~~~~~~~

Some of the dependencies of ``tqecd`` are implemented using compiled languages. This is for
example the case of the `pycryptosat <https://pypi.org/project/pycryptosat/>`_ dependency.
Pre-compiled Python packages that should be compatible with any GNU/Linux are provided
by the author, but no pre-compiled package exist for Windows or MacOS.

This means that, if you try to install the ``tqecd`` package on Windows or MacOS, a working
C++ toolchain should also be installed on your system.

Here is a list of potential issues you might encounter and how to solve them:

- `Failed building wheel for pycryptosat <https://github.com/tqec/tqec/issues/311>`_

Installation procedure
----------------------

(optional) Create a new environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is a good practice to create a specific virtual environment for the ``tqecd`` package.

One way of doing that is using the native ``venv`` package of your python installation:

.. code-block:: bash

    mkdir venvs
    python -m venv venvs/tqecd
    # On GNU/Linux and MacOS
    source venvs/tqecd/bin/activate
    # On Windows
    ## In cmd.exe
    venvs\tqecd\Scripts\activate.bat
    ## In PowerShell
    venvs\tqecd\Scripts\Activate.ps1

Install the ``tqecd`` package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``tqecd`` package is a regular Python package that can be installed using ``pip``.

It is not (yet) available on the official Python Package Index PyPI, so you will have
to manually provide the URL to install the package:

.. code-block:: bash

    python -m pip install git+https://github.com/tqec/tqecd.git

And that's it! You can test the installation by running

.. code-block:: bash

    python -c "import tqecd"

If the installation succeeded, the command should return without any message displayed.
Else, a message like

.. code-block::

    Traceback (most recent call last):
      File "<string>", line 1, in <module>
    ModuleNotFoundError: No module named 'tqecd'

should appear.
