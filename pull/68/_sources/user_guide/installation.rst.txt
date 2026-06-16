How to install ``tqecd``
========================

Requirements before installing
------------------------------

Python version
~~~~~~~~~~~~~~

The ``tqecd`` package only supports Python 3.10 and onward. If you have Python 3.9 or below,
please update your Python installation.

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
