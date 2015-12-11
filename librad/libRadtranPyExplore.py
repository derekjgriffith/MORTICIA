# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import writeLex  # this seems to load all available libradtran option definitions

# <codecell>

Options = writeLex.loadOptions()  # The files *_options.py define all the options

# <codecell>

Options # Options is a dictionary of available options

# <codecell>

dir(Options['aerosol_vulcan'])

# <codecell>

dir(Options['sza'].childs)


# <codecell>

Options['sza']['documentation']

# <codecell>

Options['sza'].tokens[0]  # Tokens are the parameters that can appear after the option keyword

# <codecell>

Options['sza'].tokens[0].datatype

# <codecell>

Options['sza'].tokens[0].get(1)

# <codecell>

writeLex.loadOptionsHelper()

# <codecell>

options['sza'].childs

# <codecell>

data, line_nos = ReadFile('UVSPEC_AEROSOL.INP')

# <codecell>

# Need another class which has the method AddValue
for line, pos in zip(data, line_nos):
    name = line[0]
    try:
        obj = Options[name]
        obj.AddValue(line[1:] if line[1:] else [])
    except KeyError:
        error = "Unknown option %s in file '%s' line %i." % \
            (repr(name), pos[0], pos[1])

# <codecell>

data

# <codecell>

line_nos

# <codecell>

basestring

# <codecell>

import widgets # Again had to edit out 

# <codecell>


