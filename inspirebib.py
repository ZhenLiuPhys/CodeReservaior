#!/usr/bin/python

# Requires MultipartPostHandler2-0.1.2 or later, which you can get from
# https://pypi.python.org/pypi/MultipartPostHandler2
# Earlier versions send too many copies of \r\n, leading to
# "internal server error".
# 
# This is Version 1.1.   Copyright 2013, 2014, 2016 Ken Olum
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import urllib
import urllib2
import MultipartPostHandler
import argparse
import os
import sys
import re

parser = argparse.ArgumentParser(description='Get bibliography information from INSPIRE')
parser.add_argument('texfile')
parser.add_argument('bibfile',nargs='?',default=None)
# Format codes: hlxu = US, hlxe = EU, hx = BibTeX
parser.add_argument('-b', dest='format', action='store_const', const='hx', help='BibTeX format')
parser.add_argument('-u', dest='format', action='store_const', const='hlxu', help='US LaTeX format')
parser.add_argument('-e', dest='format', action='store_const', const='hlxe', help='EU LaTeX format')
args = parser.parse_args()

# Default to BibTeX
if args.format == None: args.format = 'hx'

texfile = args.texfile
bibfile = args.bibfile

instream = open(texfile, "rb")

# Default bib file
if bibfile == None:
    bibfile = os.path.splitext(args.texfile)[0]+".bib"
    if os.path.exists(bibfile):
        sys.stdout.write("File "+bibfile+" exists.  Overwrite? ")
        if not raw_input().lower() in {'y', 'yes'}: exit(1)

sys.stdout.write("Contacting INSPIRE...")
sys.stdout.flush()

# First open basic page in order to obtain an access number
f = urllib2.urlopen("http://inspirehep.net/submit?ln=en&doctype=bibtex");
response = f.read()
match = re.search('<input type="hidden" name="access" value="(.*)" />',response)
if match == None: raise StandardError, "Couldn't find access number"

access =  match.expand(r'\1')

params = {'doctype': 'bibtex', 'act': 'SBI', 'OUT_FORMAT': args.format,
          'nextPg' : '', 'access': access, 'mode': 'U', 'step' : '1', 'ln': 'en',
          'curpage': '1', 'nbPg': '',
          'BibTex_input': instream}

# Make opener that knows how to do multipart encoding.
opener = urllib2.build_opener(MultipartPostHandler.MultipartPostHandler)

f = opener.open("http://inspirehep.net/submit", params)
response = f.read()
# Set DOTALL, search for multi-line preformatted string.
match = re.search(r'(?s)<pre>(.*)</pre>',response)
if match == None: raise StandardError, "Couldn't parse response from INSPIRE"
response = match.expand(r'\1')
# Now get rid of any extra <pre> tags, leaving newlines
response = re.sub('\n*<pre>\n*','\n',response)
response = re.sub('\n*</pre>\n*','\n',response)

# Get rid of warnings from INSPIRE about references it could not find
# This allows you to have a non-inspire bibliography in addition
response = re.sub(r'(?s)<div class="(nonstandard|notfound)">(.*?)</div>','\n',response)

outstream = open(bibfile,'w')
outstream.write(response)
outstream.close()

print "Bibliography written to " + bibfile
