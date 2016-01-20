MORTICIA Style Guide
=====================
 The ``MORTICIA`` style guide is a comprehensive source of all things related to ``MORTICIA`` nomenclature. This guide is developed to be compatible with python 2.7, no
 immediate consideration is given to python 3.x at this stage. ``MORTICIA`` developers suspect that any changes for python 3.x will be subtle, unobtrusive and non-detrimental
 to ``MORTICIA`` development.

To PEP-8 or not to PEP-8
========================
While PEP-8 has become the *de facto* standard for pythonic coding, it was deemed necessary by the ``MORTICIA`` developers
to make certain exceptions. 

A comprehensive PEP-8 guide can be found `here <http://legacy.python.org/dev/peps/pep-0008/#a-foolish-consistency-is-the-hobgoblin-of-little-minds>`_, while a pdf 
PEP-8 cheat sheet can be found `here <https://www.google.co.za/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwjZr-egubPKAhVH1RQKHaZKDF8QFggaMAA&url=https%3A%2F%2Fwww.pkimber.net%2Fopen%2F_downloads%2Fpep8_cheat.pdf&usg=AFQjCNGGoQ7DSCwaprDjuc356SwegyEkWA&sig2=Z-aVuunhXiWS-byXk3TGCQ>`_.

The PEP-8 guides referenced above provides a thorough encapsulation of the do's and dont's of efficient and consistent python programming, it is however not all inclusive.
There are topics which the guide does not address. This is in no way detrimental to any python development and attempting to fill these gaps are probably not an essential topic for the PEP-8 authors to address. 
Importantly, the PEP-8 is just a guide and one should not become indoctrinated with it's contents. To paraphrase `Guido Van Rossum <https://en.wikipedia.org/wiki/Guido_van_Rossum>`_ (Python's creator), it is important for ones python code 
to be consistent with PEP-8, it is more important to be consistent within a project, but it is most important to be consistent within one function or module.

The Python programming language is designed to be an *easy-to-read, easy-to-write* language. The PEP-8 style guide aims to accentuate this strong point of python by providing a means of consistency and optimal readability.

For sake of completeness the more commonly used cases of PEP-8 will be shown herein. It will be split into sub sections, after which a short section on ``MORTICIA's`` deviations will be noted. Where feasible 
a short discussion/ motivation together with an example will be provided

Code Layout
-----------

* *Indentation* : Use four space indentations. Most editors enforce this rule automatically. Spaces are preferred over tabs for indentation, tabs should be used however when adding to 
					code that already implements tabs for indentation (**consistency is key**).
				  
* *Blank lines* : Two blank lines on either side of a top level function or class definition. Method definitions inside a class are surrounded by a single blank line. Use blank lines
					sparingly inside function definitions to indicate logical sections.
				  
				  .. code-block:: python
					
					"""docstrings enclosed by triple quotes go here.
						
												
					"""
					
					import statements_go_here
					
					
					def top_level_func():
						does stuff
						return things
						
					
					class FooBar(object):
						"""Class docstrings supplied for all public classes, functions and methods
						
						"""
						
						
						def __init__(self, foo, bar):
						
							initialise things in here
							return 
							
						def internal_func(self, bar, foo):
							"""note that method this local function is surrounded by a single blank line
							"""
							perform some calculation/ data manipulation
							
							return something
					
					
					# do stuff outside the class
					
* *Line length* :	The PEP-8 style guide suggests a 79 column limit to all lines. The following is captured verbatim from the python PEP-8 style guide referenced earlier`
				  
						"The limits are chosen to avoid wrapping in editors with the window width set to 80, even if the tool places
						a marker glyph in the final column when wrapping lines. Some web based tools may not offer dynamic line 
						wrapping at all.
				  
						Some teams strongly prefer a longer line length. For code maintained exclusively or primarily by a 
						team that can reach agreement on this issue, it is okay to increase the nominal line length from 80 to 100 
						characters (effectively increasing the maximum length to 99 characters), provided that comments 
						and docstrings are still wrapped at 72 characters...."
				  
						The `MORTICIA` team have decided to conservatively follow the 99 character limit. Where possible and proves
						efficient line splitting is used.

					
						
Import Statements
------------------
This is a simple example::

    import math
    import seaborn as sns
	
Whitespace
-----------
Comments
---------
Documentation Strings (docstrings)
----------------------------------
Naming Conventions
------------------
This is an important one, the style of naming different python objects should make identifying the type of object intuitive
