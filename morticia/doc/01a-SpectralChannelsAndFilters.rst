IPython Notebook and Sphinx
===========================

To get IPython notebook, say into sphinx documentation requires converting the ipython notebook into an rst file using the following command from the command line:
ipython nbconvert hello.ipynb --to rst

This option is quite buggy at the moment. The following deficiencies have been identified:

- The code segments from the ipython notebook are converted to `code` directives. This option does not provide syntax highlighting, instead the 
  `code-block` directive needs to be manually inserted.

- In order for notebook rendered figures to be displayed, the path to the images is required to be manually altered. The `nbconvert` command creates a folder 
  in which all figures are placed, this folder must be copied to the desired location and referenced explicitly in all `image` directives.

- matplotlib objects are passed in as `parsed-literal` directives. These serve add not value to the documentation and should be deleted


Below is partially functional integration of a MORTICIA related ipython notebook. IPython notebook integration with this documentation is under development, until integration 
is stable the notebooks should be accessed from the `MORTICIA notebook repository <https://github.com/derekjgriffith/nbMORTICIA>`_  on Github 

01a-SpectralChannelsAndFilters IPython Notebook
---------------------------------------------------

.. code-block:: python

    # Perform the standard numpy and units imports
    import numpy as np
    import matplotlib.pyplot as plt
    from pint import UnitRegistry
    ureg = UnitRegistry()
    Q_ = ureg.Quantity
    import matplotlib.pyplot as plt
    %matplotlib inline

.. code-block:: python

    import sys
    sys.path.append('..')  # Remove this call if the optics module is already installed elsewhere 
    #import optics
    # This notebook is used for development/testing of the Optics module, so auto reload the Optics module if it changes
    %load_ext autoreload
    %aimport rad
    %autoreload 1

.. code-block:: python

    # Illustrate filter generation using rad.srfgen
    w,y,wn,wu = rad.srfgen(center=500, fwhm=10, shape='gauss') # Center wavelength and full width at half maximum default to nm
    plt.plot(w, y)
    plt.xlabel('Wavelength [nm]')




.. image:: 01a-SpectralChannelsAndFilters_files\01a-SpectralChannelsAndFilters_2_1.png


.. code-block:: python

    # Gaussian filter on log scale with some out-of-band leackage between 400 and 600 nm
    w,y,wn,wu = rad.srfgen(center=500, fwhm=10, shape='gauss', yedge=1e-7, oob=1e-7, wvmin=400, wvmax=600)
    plt.semilogy(wu, y)  # plot against wavelength in microns
    plt.xlabel('Wavelength [micron]')




.. parsed-literal::

    <matplotlib.text.Text at 0xb8e2518>




.. image:: 01a-SpectralChannelsAndFilters_files\01a-SpectralChannelsAndFilters_3_1.png


.. code-block:: python

    w,y,wn,wu = rad.srfgen(500, 10, shape='bartlett', wvmin=480, wvmax=520)
    plt.plot(w, y)
    plt.xlabel('Wavelength [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xbdd5b38>




.. image:: 01a-SpectralChannelsAndFilters_files\01a-SpectralChannelsAndFilters_4_1.png


.. code-block:: python

    # Bartlett (triangular) again with insertion of center flat region of width 10 nm
    w,y,wn,wu = rad.srfgen(500, 10, shape='bartlett', centerflat=10, wvmin=480, wvmax=520)
    plt.plot(w, y)
    plt.xlabel('Wavelength [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xbdad470>

.. image:: 01a-SpectralChannelsAndFilters_files\01a-SpectralChannelsAndFilters_5_1.png


.. code-block:: python

    # Generate in wavelength space (nm by default) and plot in wavenumber space
    w,y,wn,wu = rad.srfgen(500, 100, shape='welch') # center and full width at half max default to nm
    plt.plot(wn, y)  #  plot against wavenumber per cm, note how the filter is obviously skewed in wavenumber space if 
    # sufficienty wide
    plt.xlabel('Wavenumber [cm^-1]')




.. parsed-literal::

    <matplotlib.text.Text at 0xc2e0a58>




.. image:: 01a-SpectralChannelsAndFilters_files\01a-SpectralChannelsAndFilters_6_1.png


.. code-block:: python

    # Similar, but this time, generate the filter in wavenumber scale and also plot in wavenumber scale
    w,y,wn,wu = rad.srfgen(20000.0, 2000.0, shape='welch', units='cm^-1')
    plt.plot(wn, y)  #  plot against wavenumber per cm
    plt.xlabel('Wavenumber [per cm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xc48b2b0>




.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_7_1.png


.. code-block:: python

    w,y,wn,wu = rad.srfgen(500, 10, shape='cosine')
    plt.plot(w, y)
    plt.xlabel('Wavenumber [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xc5bfe10>




.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_8_1.png


.. code-block:: python

    w,y,wn,wu = rad.srfgen(500, 10, shape='cos^2', wvmin=480, wvmax=520, oob=0.001)
    plt.plot(w, y)
    plt.xlabel('Wavelength [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xaf404a8>




.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_9_1.png


.. code-block:: python

    # Same, plotted on a log scale
    plt.semilogy(w, y)
    plt.xlabel('Wavelength [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xcbc0b70>




.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_10_1.png


.. code-block:: python

    # Same as previous, but open a central flat region of 10 nm width
    # Note that the full width at half max is now 20 nm (fwhm + centerflat)
    w,y,wn,wu = rad.srfgen(500, 10, shape='cos^2', centerflat=10, wvmin=480, wvmax=520, oob=0.001)
    plt.semilogy(w, y)
    plt.xlabel('Wavelength [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xd299d30>




.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_11_1.png


.. code-block:: python

    w,y,wn,wu = rad.srfgen(605, 10, shape='tophat', wvmin=580, wvmax=640)
    plt.plot(w, y)
    plt.xlabel('Wavelength [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xd58ada0>




.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_12_1.png


.. code-block:: python

    # Opening a centre flat region also works with tophats (box) although you could just increase fwhm
    w,y,wn,wu = rad.srfgen(605, 10, shape='tophat', centerflat=20, wvmin=580, wvmax=640)
    plt.plot(w, y)
    plt.xlabel('Wavelength [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xd5eeac8>




.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_13_1.png


.. code-block:: python

    # Show the tophat function, which defines a tophat using only a few points
    w,y,wn,wu = rad.tophat(605, 10, delta=0.001, wvmin=350, wvmax=750)
    plt.plot(w,y)
    plt.xlabel('Wavelength [nm]')




.. parsed-literal::

    <matplotlib.text.Text at 0xd782e10>




.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_14_1.png


.. code-block:: python

    # Create 3 MODTRAN-style flt filters/SRFs with different postions, widths and shapes
    filt = rad.Flt('My Special Filters', filterheaders = ['a', 'b', 'c'], centers = [500, 600, 700], fwhms = [10, 20, 30],
                   shapes=['gauss', 'cos^2', 'welch'])

.. code-block:: python

    # Write the filters/SRFs in MODTRAN .flt format to a text file
    filt.write('SpecialFilters')

.. code-block:: python

    # Plot the filters
    filt.plot()



.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_17_0.png


.. code-block:: python

    # Create some filters with just a few points
    filt2 = rad.Flt('Hand Filters', filterheaders=['a', 'b'], filters=[np.array([[300, 0.5],[400, 1.0]]),
                                                                       np.array([[300, 0.3],[400, 0.2]])])

.. code-block:: python

    filt2.plot()



.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_19_0.png


.. code-block:: python

    filt2




.. parsed-literal::

    N Hand Filters
    a
    3.000000000000000000e+02 5.000000000000000000e-01 3.333333333333333576e+04
    4.000000000000000000e+02 1.000000000000000000e+00 2.500000000000000000e+04
    b
    3.000000000000000000e+02 2.999999999999999889e-01 3.333333333333333576e+04
    4.000000000000000000e+02 2.000000000000000111e-01 2.500000000000000000e+04



.. code-block:: python

    CIE = rad.Flt('CIEXYZ2.flt')

.. code:: python

    CIE.plot()



.. image:: 01a-SpectralChannelsAndFilters_files%5C01a-SpectralChannelsAndFilters_22_0.png


.. code-block:: python

    CIE




.. parsed-literal::

    N CIE 1931 Color Observer Functions, 2 deg
    X (2 deg)
      380.000000   0.001400   26315.789474
      385.000000   0.002200   25974.025974
      .
	  .
	  .
	  735.000000   0.000000   13605.442177
      740.000000   0.000000   13513.513514
      745.000000   0.000000   13422.818792
      750.000000   0.000000   13333.333333
      755.000000   0.000000   13245.033113
      760.000000   0.000000   13157.894737
      765.000000   0.000000   13071.895425
      770.000000   0.000000   12987.012987
      775.000000   0.000000   12903.225806
      780.000000   0.000000   12820.512821



.. code-block:: python

    Aviris = rad.Flt('')





.. parsed-literal::

    Nanometer data for AVIRIS sensor (assumes Gaussian with maximum response of 1)
    CENTER:  373.40 NM   FWHM:  9.90 NM
      357.773500   0.001000   27950.644751
      2521.181200   0.001530   3966.394799
      2521.759300   0.001000   3965.485524



.. code-block:: python

    Aviris.write('NewAvirisData.flt')

.. code-block:: python

    import numpy as np
    from astropy.nddata import NDData

.. code-block:: python

    a = np.array([1.1,2.2,3.3])

.. code-block:: python

    b = NDData(a)

.. code-block:: python

    c = NDData([5.5,6.6,7.7])

.. code-block:: python

    c




.. parsed-literal::

    NDData([ 5.5,  6.6,  7.7])



.. code-block:: python

    b




.. parsed-literal::

    NDData([ 1.1,  2.2,  3.3])



.. code-block:: python

    b+c


::


    ---------------------------------------------------------------------------
    TypeError                                 Traceback (most recent call last)

    <ipython-input-16-9533ace22465> in <module>()
    ----> 1 b+c
    

    TypeError: unsupported operand type(s) for +: 'NDData' and 'NDData'


.. code-block:: python

    import numpy as np
    import pandas as pd
    import xray


::


    ---------------------------------------------------------------------------
    ImportError                               Traceback (most recent call last)

    <ipython-input-2-36929e266e81> in <module>()
          1 import numpy as np
    ----> 2 import pandas as pd
          3 import xray
    

    C:\Anaconda\lib\site-packages\pandas\__init__.py in <module>()
         11                       "pandas from the source directory, you may need to run "
         12                       "'python setup.py build_ext --inplace' to build the C "
    ---> 13                       "extensions first.".format(module))
         14 
         15 from datetime import datetime
    

    ImportError: C extension: hashtable not built. If you want to import pandas from the source directory, you may need to run 'python setup.py build_ext --inplace' to build the C extensions first.


