Installation
============

Installing `mofstructure` is a straightforward process and you have two options depending on your preferred method: using `pip` or cloning the repository directly from GitHub. Below are the step-by-step instructions for both methods.

Option 1: Installing via pip
----------------------------

The easiest and most convenient way to install `mofstructure` is by using `pip`. This method ensures that latest stable release of the module along with any dependencies are installed.

To install `mofstructure` using `pip`, simply open your terminal or command prompt and execute the following command:

.. code-block:: bash

   pip install mofstructure

This command will automatically download and install the `mofstructure` package from the Python Package Index (PyPI). After the installation is complete, you will be able to import and use the `mofstructure` module in your Python projects.

Option 2: Installing via GitHub Repository
------------------------------------------

If you prefer to work with the latest development version of `mofstructure` or contribute to its development or use a different python version, you can install the module directly from its GitHub repository. This option allows you to access the most up-to-date code, including any recent changes or experimental features that may not yet be available through `pip`.

Follow the steps below to clone the repository and install `mofstructure` locally:

1. **Clone the Repository:**

   First, clone the `mofstructure` repository from GitHub to your local machine using the following command:

   .. code-block:: bash

      git clone https://github.com/bafgreat/mofstructure.git mofstructure

   This command creates a new directory named `mofstructure` in your current working directory, containing all the files from the repository.

2. **Navigate to the Directory:**

   Change your working directory to the newly cloned `mofstructure` folder:

   .. code-block:: bash

      cd mofstructure

   This step ensures that you are in the correct directory where the `setup.py` file is located, which is necessary for the installation process.

3. **Install the Module:**

   Finally, install `mofstructure` by running the following command:

   .. code-block:: bash

      pip install .

   The `pip install .` command tells `pip` to install the package from the current directory, which contains the latest version of the module's source code.

Once the installation is complete, you will have `mofstructure` available in your Python environment, and you can begin using it. This method is particularly useful if you plan to make modifications to the code or want to stay up-to-date with the latest commits.
