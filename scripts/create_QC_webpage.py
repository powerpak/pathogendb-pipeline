#!/usr/bin/env python

import subprocess
import sys
import numpy
from collections import namedtuple as ntup
from optparse import OptionParser
import os
import datetime
import networkx as nx
import shlex
import shutil


# takes RGB values and returns a colour string
def colorstr(rgb):
    return "#%02x%02x%02x" % (rgb[0],rgb[1], rgb[2])

# takes hue, saturation and lightness values and returns a color string
def hsl_to_colorstr(h, s, l):
    c = (1 - abs(2*l - 1)) * s
    x = c * (1 - abs(h *1.0 / 60 % 2 - 1))
    m = l - c/2
    if h < 60:
        r, g, b = c + m, x + m, 0 + m
    elif h < 120:
        r, g, b = x + m, c+ m, 0 + m
    elif h < 180:
        r, g, b = 0 + m, c + m, x + m
    elif h < 240:
        r, g, b, = 0 + m, x + m, c + m
    elif h < 300:
        r, g, b, = x + m, 0 + m, c + m
    else:
        r, g, b, = c + m, 0 + m, x + m
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)
    return "#%02x%02x%02x" % (r,g,b)

# This script needs this CSS template to create the website -- just put it here for safekeeping
css_template = '''body{
margin: 0;
padding: 0;
font: normal 90% "Lucida Grande", "Trebuchet MS", Verdana, Helvetica, sans-serif;
}

#framecontent{
position: absolute;
top: 0;
bottom: 0;
left: 0;
width: 200px; /*Width of frame div*/
height: 100%;
overflow: auto; /*Disable scrollbars. Set to "scroll" to enable*/
background: #2fa4e7;
color: white;
font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
}

h1 {
text-align: center;
font-weight: bold;
}

#framecontent a:link {
    color: #ffffcc;
}
#framecontent a:visited {
    color: #ffffcc;
}
#framecontent a:hover {
    color: #cc0066;
}




#maincontent{
position: fixed;
top: 0;
left: 200px; /*Set left value to WidthOfFrameDiv*/
right: 0;
bottom: 0;
overflow: auto;
background: #fff;
}


#maincontent h1{
color: #157ab5;
}

.innertube{
margin: 15px; /*Margins for inner DIV inside each DIV (to provide padding)*/
}

* html body{ /*IE6 hack*/
padding: 0 0 0 300px; /*Set value to (0 0 0 WidthOfFrameDiv)*/
}

* html #maincontent{ /*IE6 hack*/
height: 100%;
width: 100%;
}


table {
	width: 100%;
	border-collapse: collapse;
        font-size: 12px;
}

table a:hover {
  color: #157ab5;
  text-decoration: underline;
}
tr {
        background: #f9f9f9;
        color: #333333
}
tr:nth-of-type(odd) {
	background: #f5f5f5;
}
th {
	background: #d4d4d4;
	color: #157ab5;
	font-weight: bold;
}
td, th {
	padding: 6px;
	border: 1px solid #ccc;
	text-align: left;
}
thead {
    display: inline-block;
    height: 20px;
}
tbody {
    height: 200px;
    display: inline-block;
    overflow: auto;'''


# class for creating scalable vector graphics
class scalableVectorGraphics:

    def __init__(self, height, width):
        self.height = height
        self.width = width
        self.out = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   height="%d"
   width="%d"
   id="svg2"
   version="1.1"
   inkscape:version="0.48.4 r9939"
   sodipodi:docname="easyfig">
  <metadata
     id="metadata122">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title>Easyfig</dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <defs
     id="defs120" />
  <sodipodi:namedview
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1"
     objecttolerance="10"
     gridtolerance="10"
     guidetolerance="10"
     inkscape:pageopacity="0"
     inkscape:pageshadow="2"
     inkscape:window-width="640"
     inkscape:window-height="480"
     id="namedview118"
     showgrid="false"
     inkscape:zoom="0.0584"
     inkscape:cx="2500"
     inkscape:cy="75.5"
     inkscape:window-x="55"
     inkscape:window-y="34"
     inkscape:window-maximized="0"
     inkscape:current-layer="svg2" />
  <title
     id="title4">Easyfig</title>
  <g
     style="fill-opacity:1.0; stroke:black; stroke-width:1;"
     id="g6">''' % (self.height, self.width)

    def drawLine(self, x1, y1, x2, y2, th=1, cl=(0, 0, 0), alpha = 1.0):
        self.out += '  <line x1="%d" y1="%d" x2="%d" y2="%d"\n        stroke-width="%d" stroke="%s" stroke-opacity="%f" stroke-linecap="butt" />\n' % (x1, y1, x2, y2, th, colorstr(cl), alpha)

    def writesvg(self, filename):
        outfile = open(filename, 'w')
        outfile.write(self.out + ' </g>\n</svg>')
        outfile.close()

    def drawRightArrow(self, x, y, wid, ht, fc, oc=(0,0,0), lt=1):
        if lt > ht /2:
            lt = ht/2
        x1 = x + wid
        y1 = y + ht/2
        x2 = x + wid - ht / 2
        ht -= 1
        if wid > ht/2:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x, y+ht/4, x2, y+ht/4,
                                                                                                x2, y, x1, y1, x2, y+ht,
                                                                                                x2, y+3*ht/4, x, y+3*ht/4)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x, y, x, y+ht, x + wid, y1)

    def drawLeftArrow(self, x, y, wid, ht, fc, oc=(0,0,0), lt=1):
        if lt > ht /2:
            lt = ht/2
        x1 = x + wid
        y1 = y + ht/2
        x2 = x + ht / 2
        ht -= 1
        if wid > ht/2:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x1, y+ht/4, x2, y+ht/4,
                                                                                                x2, y, x, y1, x2, y+ht,
                                                                                                x2, y+3*ht/4, x1, y+3*ht/4)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x, y1, x1, y+ht, x1, y)

    def drawBlastHit(self, x1, y1, x2, y2, x3, y3, x4, y4, fill=(0, 0, 255), lt=2, alpha=0.1):
        self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d" fill-opacity="%f"\n' % (colorstr(fill), colorstr((0,0,0)), lt, alpha)
        self.out += '           points="%d,%d %d,%d %d,%d %d,%d" />\n' % (x1, y1, x2, y2, x3, y3, x4, y4)

    def drawGradient(self, x1, y1, wid, hei, minc, maxc):
        self.out += '  <defs>\n    <linearGradient id="MyGradient" x1="0%" y1="0%" x2="0%" y2="100%">\n'
        self.out += '      <stop offset="0%%" stop-color="%s" />\n' % colorstr(maxc)
        self.out += '      <stop offset="100%%" stop-color="%s" />\n' % colorstr(minc)
        self.out += '    </linearGradient>\n  </defs>\n'
        self.out += '  <rect fill="url(#MyGradient)" stroke-width="0"\n'
        self.out += '        x="%d" y="%d" width="%d" height="%d"/>\n' % (x1, y1, wid, hei)

    def drawGradient2(self, x1, y1, wid, hei, minc, maxc):
        self.out += '  <defs>\n    <linearGradient id="MyGradient2" x1="0%" y1="0%" x2="0%" y2="100%">\n'
        self.out += '      <stop offset="0%%" stop-color="%s" />\n' % colorstr(maxc)
        self.out += '      <stop offset="100%%" stop-color="%s" />\n' % colorstr(minc)
        self.out += '    </linearGradient>\n</defs>\n'
        self.out += '  <rect fill="url(#MyGradient2)" stroke-width="0"\n'
        self.out += '        x="%d" y="%d" width="%d" height="%d" />\n' % (x1, y1, wid, hei)

    def drawOutRect(self, x1, y1, wid, hei, fill, lt=1, alpha=0.3):
        self.out += '  <rect fill="%s" stroke-width="%d" fill-opacity="%f"\n' % (colorstr(fill), lt, alpha)
        self.out += '        x="%d" y="%d" width="%d" height="%d" />\n' % (x1, y1, wid, hei)

    def drawRightFrame(self, x, y, wid, ht, lt, frame, fill):
        if lt > ht /2:
            lt = ht /2
        if frame == 1:
            y1 = y + ht/2
            y2 = y + ht * 3/8
            y3 = y + ht * 1/4
        elif frame == 2:
            y1 = y + ht * 3/8
            y2 = y + ht * 1/4
            y3 = y + ht * 1/8
        elif frame == 0:
            y1 = y + ht * 1/4
            y2 = y + ht * 1/8
            y3 = y + 1
        x1 = x
        x2 = x + wid - ht/8
        x3 = x + wid
        if wid > ht/8:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x1, y1, x2, y1, x3, y2, x2, y3, x1, y3)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x1, y1, x3, y2, x1, y3)

    def drawRightFrameRect(self, x, y, wid, ht, lt, frame, fill):
        if lt > ht /2:
            lt = ht /2
        if frame == 1:
            y1 = y + ht / 4
        elif frame == 2:
            y1 = y + ht /8
        elif frame == 0:
            y1 = y + 1
        hei = ht /4
        x1 = x
        self.out += '  <rect fill="%s" stroke-width="%d"\n' % (colorstr(fill), lt)
        self.out += '        x="%d" y="%d" width="%d" height="%d" />\n' % (x1, y1, wid, hei)

    def drawLeftFrame(self, x, y, wid, ht, lt, frame, fill):
        if lt > ht /2:
            lt = ht /2
        if frame == 1:
            y1 = y + ht
            y2 = y + ht * 7/8
            y3 = y + ht * 3/4
        elif frame == 2:
            y1 = y + ht * 7/8
            y2 = y + ht * 3/4
            y3 = y + ht * 5/8
        elif frame == 0:
            y1 = y + ht * 3/4
            y2 = y + ht * 5/8
            y3 = y + ht / 2
        x1 = x + wid
        x2 = x + ht/8
        x3 = x
        if wid > ht/8:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
            self.out += '           points="%d,%d %d,%d %d,%d %d,%d %d,%d" />\n' % (x1, y1, x2, y1, x3, y2, x2, y3, x1, y3)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
            self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x1, y1, x3, y2, x1, y3)

    def drawLeftFrameRect(self, x, y, wid, ht, lt, frame, fill):
        if lt > ht /2:
            lt = ht /2
        if frame == 1:
            y1 = y + ht * 3/4
        elif frame == 2:
            y1 = y + ht * 5/8
        elif frame == 0:
            y1 = y + ht / 2
        hei = ht /4
        x1 = x
        self.out += '  <rect fill="%s" stroke-width="%d"\n' % (colorstr(fill), lt)
        self.out += '        x="%d" y="%d" width="%d" height="%d" />\n' % (x1, y1, wid, hei)

    def drawPointer(self, x, y, ht, lt, fill):
        x1 = x - int(round(0.577350269 * ht/2))
        x2 = x + int(round(0.577350269 * ht/2))
        y1 = y + ht/2
        y2 = y + 1
        self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%d"\n' % (colorstr(fill), colorstr((0, 0, 0)), lt)
        self.out += '           points="%d,%d %d,%d %d,%d" />\n' % (x1, y2, x2, y2, x, y1)

    def drawDash(self, x1, y1, x2, y2, exont):
        self.out += '  <line x1="%d" y1="%d" x2="%d" y2="%d"\n' % (x1, y1, x2, y2)
        self.out += '       style="stroke-dasharray: 5, 3, 9, 3"\n'
        self.out += '       stroke="#000" stroke-width="%d" />\n' % exont

    def drawPolygon(self, x_coords, y_coords, colour=(0,0,255)):
        self.out += '  <polygon points="'
        for i,j in zip(x_coords, y_coords):
            self.out += str(i) + ',' + str(j) + ' '
        self.out += '"\nstyle="fill:%s;stroke=none" />\n'  % colorstr(colour)
    def writeString(self, thestring, x, y, size, ital=False, bold=False, rotate=0, justify='left'):
        if rotate != 0:
            x, y = y, x
        self.out += '  <text\n'
        self.out += '    style="font-size:%dpx;font-style:normal;font-weight:normal;z-index:10\
;line-height:125%%;letter-spacing:0px;word-spacing:0px;fill:#999999;fill-opacity:1;stroke:none;font-family:Sans"\n' % size
        if justify == 'right':
            self.out += '    text-anchor="end"\n'
        if rotate == 1:
            self.out += '    x="-%d"\n' % x
        else:
            self.out += '    x="%d"\n' % x
        if rotate == -1:
            self.out += '    y="-%d"\n' % y
        else:
            self.out += '    y="%d"\n' % y
        self.out += '    sodipodi:linespacing="125%"'
        if rotate == -1:
            self.out += '\n    transform="matrix(0,1,-1,0,0,0)"'
        if rotate == 1:
            self.out += '\n    transform="matrix(0,-1,1,0,0,0)"'
        self.out += '><tspan\n      sodipodi:role="line"\n'
        if rotate == 1:
            self.out += '      x="-%d"\n' % x
        else:
            self.out += '      x="%d"\n' % x
        if rotate == -1:
            self.out += '      y="-%d"' % y
        else:
            self.out += '      y="%d"' % y
        if ital and bold:
            self.out += '\nstyle="font-style:italic;font-weight:bold"'
        elif ital:
            self.out += '\nstyle="font-style:italic"'
        elif bold:
            self.out += '\nstyle="font-style:normal;font-weight:bold"'
        self.out += '>' + thestring + '</tspan></text>\n'

# class for reading sam files - emulates pysam - except runs on windows/doesn't need a library
class readSam:
    def __init__(self, sam_file):
        self.header = ''
        self.references = []
        self.lengths = []
        self.read = ntup('blast', 'pos pnext rname rnext is_reverse mate_is_reverse is_read1 is_unmapped mate_is_unmapped line cigar align_end read_length name large_inserts large_dels start_is_trimmed end_is_trimmed secondary_alignment')
        self.sam = open(sam_file)
        line = self.sam.readline()
        lastline = None
        while line.startswith('@'):
            self.header += line
            if line.startswith('@SQ'):
                for i in line.split():
                    if i.startswith('SN:'):
                        self.references.append(i[3:])
                    elif i.startswith('LN:'):
                        self.lengths.append(int(i[3:]))
            lastline = line
            line = self.sam.readline()
        self.sam.seek(0)
        if not lastline is None:
            getit = True
            while getit:
                line = self.sam.readline()
                if line == lastline:
                    getit = False

    def __iter__(self):
        return self

    def next(self):
        line = self.sam.readline()
        if line == '':
            raise StopIteration
        name, flag, rname, pos, mapq, cigar, rnext, pnext, tlength, sequence = line.split()[:10]
        if rnext == '=' or rnext == '*':
            rnext = rname
        flag = bin(int(flag)).zfill(12)
        currdigit = ''
        readlength = 0
        alignlength = 0
        inserts = []
        deletions = []
        trimmed_start = None
        for i in cigar:
            if i == '*':
                pass
            elif i.isdigit():
                currdigit += i
            else:
                thenum = int(currdigit)
                currdigit = ''
                if i == 'M' or i == '=' or i == 'X':
                    readlength += thenum
                    alignlength += thenum
                elif i == 'I':
                    readlength += thenum
                    if thenum > large_gap_min:
                        inserts.append((int(pos) + alignlength - 1, thenum))
                elif i == 'D':
                    alignlength += thenum
                    if thenum > large_gap_min:
                        deletions.append((int(pos) + alignlength - 1, thenum))
                elif i == 'H' or i == 'S':
                    readlength += thenum
                    if trimmed_start is None:
                        trimmed_start = True
                elif i == 'N':
                    alignlength += thenum
                    if thenum > large_gap_min:
                        deletions.append((int(pos) + alignlength - 1, thenum))
                if trimmed_start is None:
                    trimmed_start = False
        if cigar[-1] == 'H' or cigar[-1] == 'S':
            trimmed_end = True
        else:
            trimmed_end = False
        align_end = int(pos) + alignlength - 1
        read = self.read(int(pos), int(pnext), rname, rnext, flag[-5] == '1', flag[-6] == '1', flag[-7] == '1', flag[-3] == '1', flag[-4] == '1', line, cigar, align_end, readlength, name, inserts, deletions, trimmed_start, trimmed_end, flag[-9] == '1')
        return read



# create the graph webpages - also create a bigwig file of the graphs
def draw_graph(options, header, footer):
    sam_filename = options.output_folder + '/alignment.sam'
    sam = readSam(sam_filename)
    bin_size = options.bin_size
    bin_step = options.bin_step
    flag_ratio = options.flag_ratio
    out_directory = options.output_folder
    out_cov = {}
    out_flag = {}
    chrom_size = {}
    for refnum, reference in enumerate(sam.references):
        html_name = reference[:6]
        sam = readSam(sam_filename)
        bin_num = sam.lengths[refnum] / bin_step + 1
        forward_start = numpy.zeros(bin_num, dtype=numpy.int)
        forward_end = numpy.zeros(bin_num, dtype=numpy.int)
        forward_start_clipped = numpy.zeros(bin_num, dtype=numpy.int)
        forward_end_clipped = numpy.zeros(bin_num, dtype=numpy.int)
        forward_through = numpy.zeros(bin_num, dtype=numpy.int)
        reverse_start = numpy.zeros(bin_num, dtype=numpy.int)
        reverse_end = numpy.zeros(bin_num, dtype=numpy.int)
        reverse_start_clipped = numpy.zeros(bin_num, dtype=numpy.int)
        reverse_end_clipped = numpy.zeros(bin_num, dtype=numpy.int)
        reverse_through = numpy.zeros(bin_num, dtype=numpy.int)
        large_deletions = numpy.zeros(bin_num, dtype=numpy.int)
        large_insertions = numpy.zeros(bin_num, dtype=numpy.int)
        coverage_array = numpy.zeros(bin_num)
        if sam.lengths[refnum] <= 300000:
            x_axis = numpy.arange(-sam.lengths[refnum]/3, sam.lengths[refnum]/3*2 + 1, bin_step)
            chrom_size[reference] = sam.lengths[refnum]/3
        else:
            x_axis = numpy.arange(-100000, sam.lengths[refnum] -100000 + 1, bin_step)
            chrom_size[reference] = sam.lengths[refnum] - 200000
        last_read, last_cigar = None, None
        for read in sam:
            if read.rname == reference:
                gotten_indels = set()
                for i in read.large_dels:
                    for j in range(max([0, (i[0] - bin_size)/bin_step]), i[0]/bin_step + 1):
                        if not j in gotten_indels:
                            large_deletions[j] += 1
                            gotten_indels.add(j)
                gotten_indels = set()
                for i in read.large_inserts:
                    for j in range(max([0, (i[0] - bin_size)/bin_step]), i[0]/bin_step + 1):
                        if not j in gotten_indels:
                            large_insertions[j] += 1
                            gotten_indels.add(j)
                start = read.pos
                end = read.align_end
                for i in range(max([0, (start - bin_size)/bin_step]), end/bin_step + 1):
                    if not read.secondary_alignment or (read.name == last_read and read.cigar == last_cigar):
                        if i * bin_step <= start and i * bin_step + bin_size > start and i * bin_step + bin_size <= end:
                            if read.start_is_trimmed and read.is_reverse:
                                reverse_start_clipped[i] += 1
                            elif read.start_is_trimmed:
                                forward_start_clipped[i] += 1
                            elif read.is_reverse:
                                reverse_start[i] += 1
                            else:
                                forward_start[i] += 1
                            coverage_array[i] += ((i * bin_step + bin_size) - start) * 1.0 / bin_size
                        if start < i * bin_step and end > i * bin_step + bin_size:
                            if read.is_reverse:
                                reverse_through[i] += 1
                            else:
                                forward_through[i] += 1
                            coverage_array[i] += 1
                        if i * bin_step <= end and i * bin_step + bin_size > end and start < i * bin_step:
                            if read.end_is_trimmed and read.is_reverse:
                                reverse_end_clipped[i] += 1
                            elif read.end_is_trimmed:
                                forward_end_clipped[i] += 1
                            elif read.is_reverse:
                                reverse_end[i] += 1
                            else:
                                forward_end[i] += 1
                            coverage_array[i] += (end - i * bin_step) * 1.0 / bin_size
                        if i * bin_step <= start and i * bin_step + bin_size > end:
                            coverage_array[i] += (end - start) * 1.0 / bin_size
                        last_read = read.name
                        last_cigar = read.cigar
            else:
                if not read.secondary_alignment:
                    last_read = read.name
                    last_cigar = read.cigar
        flag_end_clipped = -10
        flag_end_unclipped = -10
        flag_end_indel = -10
        clipped_flag = []
        unclipped_flag = []
        indel_flag = []
        if sam.lengths[refnum] >= 300000:
            offset = 100000
        else:
            offset = sam.lengths[refnum] / 3
        for ft,rt,fs,rs,fe,re,fsc,rsc,fec,rec,li,ld,cov,count in zip(forward_through, reverse_through, forward_start, reverse_start, forward_end, reverse_end, forward_start_clipped,
                          reverse_start_clipped, forward_end_clipped, reverse_end_clipped, large_insertions, large_deletions, coverage_array, range(len(forward_through))):
             if fec >= flag_ratio * (ft + fec) or rec >= flag_ratio * (rt + rec) or fsc >= flag_ratio * (ft + fsc) or rsc >= flag_ratio * (rsc + rt):
                 if flag_end_clipped + 1 == count:
                     flag_end_clipped = count
                 else:
                     if not flag_end_clipped == -10:
                         clipped_flag.append((flag_start_clipped * bin_step - offset, flag_end_clipped * bin_step - offset + bin_size,
                                              'Greater than ' + str(int(flag_ratio * 100)) + ' % of reads Trimmed in bin.'))
                     flag_start_clipped = count
                     flag_end_clipped = count
             if fe >= flag_ratio * (ft + fe) or re >= flag_ratio * (rt + re) or fs >= flag_ratio * (fs + ft) or rs >= flag_ratio * (rs + rt):
                 if flag_end_unclipped + 1 == count:
                     flag_end_unclipped = count
                 else:
                     if not flag_end_unclipped == -10:
                         unclipped_flag.append((flag_start_unclipped * bin_step - offset, flag_end_unclipped * bin_step - offset + bin_size,
                                                'Greater than ' + str(int(flag_ratio * 100)) + ' % of reads end (untrimmed) in bin.'))
                     flag_start_unclipped = count
                     flag_end_unclipped = count
             if li >= flag_ratio * cov or ld >= flag_ratio * cov:
                 if flag_end_indel + 1 == count:
                     flag_end_indel = count
                 else:
                     if not flag_end_indel == -10:
                         indel_flag.append((flag_start_indel * bin_step - offset, flag_end_indel * bin_step - offset + bin_size,
                                            'Greater than ' + str(int(flag_ratio * 100)) + ' % of reads have an indel > ' + str(large_gap_min)))
                     flag_start_indel = count
                     flag_end_indel = count
        html_out = open(out_directory + '/qc_website/graphs/' + str(html_name) + '_graphs.html', 'w')
        html_out.write(header)
        html_out.write('  <script type="text/javascript">\n'
                       '  window.onload = function () {\n'
                       '  var chart1 = new CanvasJS.Chart("chartContainer1",\n'
                       '    {\n'
                       '      zoomEnabled: true,\n'
                       '      title:{\n'
                       '      text: "Total number of each read type in bin - bin size: ' + str(bin_size) + ', bin step: ' + str(bin_step) + ')",\n'
                       '      fontSize: 24\n'
                       '      },\n'
                       '      axisX: {\n'
                       '      title: "Position in genome",\n'
                       '      titleFontSize: 16,\n')
        if clipped_flag != [] or unclipped_flag != []:
            html_out.write('      stripLines:[\n')
            for values in clipped_flag[:-1]:
                html_out.write('      {\n'
                               '        startValue: ' + str(values[0]) + ',\n'
                               '        endValue: ' + str(values[1]) + ',\n'
                               '        color: "#16CBEF"\n'
                               '      },\n')
            if clipped_flag != []:
                html_out.write('      {\n'
                               '        startValue: ' + str(clipped_flag[-1][0]) + ',\n'
                               '        endValue: ' + str(clipped_flag[-1][1]) + ',\n'
                               '        color: "#16CBEF"\n'
                               '      }')
                if unclipped_flag != []:
                    html_out.write(',')
                else:
                    html_out.write('\n')
            for values in unclipped_flag[:-1]:
                html_out.write('      {\n'
                               '        startValue: ' + str(values[0]) + ',\n'
                               '        endValue: ' + str(values[1]) + ',\n'
                               '        color: "#53A01F"\n'
                               '      },\n')
            if unclipped_flag != []:
                html_out.write('      {\n'
                               '        startValue: ' + str(unclipped_flag[-1][0]) + ',\n'
                               '        endValue: ' + str(unclipped_flag[-1][1]) + ',\n'
                               '        color: "#53A01F"\n'
                               '      }\n')
            html_out.write('       ]\n  ')

        html_out.write('      },\n'
                       '      axisY: {\n'
                       '      title: "Number of reads",\n'
                       '      titleFontSize: 16\n'
                       '      },\n'
                       '       data: [')
        for the_data, leg_lab, colour in zip([forward_through, reverse_through, forward_start, reverse_start, forward_end, reverse_end, forward_start_clipped,
                         reverse_start_clipped, forward_end_clipped, reverse_end_clipped], ['Read spans bin (Forward)',
                          'Read spans bin (Reverse)', 'Read starts in bin (F)', 'Read starts in bin (R)',
                          'Read terminates in bin (F)', 'Read terminates in bin (R)', 'Read start clipped in bin (F)',
                          'Read start clipped in bin (R)', 'Read end clipped in bin (F)', 'Read end clipped in bin (R)'],
                          ["#BE4200", "#EE16B8", "#38081F", "#4282F6", "#45FAA5", "#763181",
                              "#234D65", "#F1F272", "#D15AF4", "#CBC22F"]):
            html_out.write('\n{\n'
                           '         type: "line",\n'
                           '         showInLegend: true,\n'
                           '         legendText: "' + leg_lab + '",\n'
                           '         color: "' + colour + '",\n'
                           '         markerType: "none",\n'
                           '         dataPoints: [\n')
            if refnum == 0:
                bg_out = open(options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() + '.wig', 'w')
                bgt_out = open(options.output_folder + '/bigwig/' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() + '.bwt', 'w')
                bgt_out.write('track type=bigwig bigDataUrl=https://vanbah01.u.hpc.mssm.edu/igb/' + options.assembly_name +'/'
                              + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower()
                              + '.bw name=' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() +
                              'color=0,0,200 altColor=0,200,0 autoScale=on alwaysZero=on graphType=bar yLineMark=10 yLineOnOff=on\n')
                bgt_out.close()
            else:
                bg_out = open(options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() + '.wig', 'a')
            bg_out.write('fixedStep chrom=' + reference + ' start=1 step=' + str(bin_step) + ' span=' + str(bin_step) + '\n')
            for value in range(0, len(x_axis) - 1):
                html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')
                if x_axis[value] >= 0 and x_axis[value] < chrom_size[reference] - bin_step:
                    bg_out.write(str(the_data[value]) + '\n')
            bg_out.close()
            html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
            html_out.write('        ]\n      }')
            if not the_data is reverse_end_clipped:
                html_out.write(',\n')
        html_out.write('      ],\n'
                       '      rangeChanged: syncHandler\n'
                       '    });\n'
                       '    chart1.render();\n'
                       '  var chart2 = new CanvasJS.Chart("chartContainer2",\n'
                       '    {\n'
                       '      zoomEnabled: true,\n'
                       '      title:{\n'
                       '      text: "Proportion of each read type in bin - bin size: ' + str(bin_size) + ', bin step: ' + str(bin_step) + ')",\n'
                       '      fontSize: 24\n'
                       '      },\n'
                       '      axisX: {\n'
                       '      title: "Position in genome",\n'
                       '      titleFontSize: 16\n'
                       '},\n'
                       '      axisY: {\n'
                       '      title: "Number of reads",\n'
                       '      titleFontSize: 16\n'
                       '},\n'
                       '       data: [')
        for the_data, leg_lab in zip([forward_through, reverse_through, forward_start, reverse_start, forward_end, reverse_end, forward_start_clipped,
                         reverse_start_clipped, forward_end_clipped, reverse_end_clipped], ['Read spans bin (Forward)',
                          'Read spans bin (Reverse)', 'Read starts in bin (F)', 'Read starts in bin (R)',
                          'Read terminates in bin (F)', 'Read terminates in bin (R)', 'Read start clipped in bin (F)',
                          'Read start* clipped in bin (R)', 'Read end clipped in bin (F)', 'Read end clipped in bin (R)']):
            html_out.write('         {\n'
                           '         type: "stackedArea100",\n'
                           '         showInLegend: true,\n'
                           '         legendText: "' + leg_lab + '",\n'
                           '         markerType: "none",\n'
                           '         legendMarkerType: "square",\n'
                           '        dataPoints: [\n')
            for value in range(0, len(x_axis) - 1):
                html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')

            html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
            html_out.write('        ]\n      }')
            if not the_data is reverse_end_clipped:
                html_out.write(',')
        html_out.write('      ],\n'
                       '      rangeChanged: syncHandler\n'
                       '   });\n'
                       '    chart2.render();\n'
                       '  var chart3 = new CanvasJS.Chart("chartContainer3",\n'
                       '    {\n'
                       '      zoomEnabled: true,\n'
                       '      title:{\n'
                       '      text: "Number of large indels in each bin - bin size: ' + str(bin_size) + ', bin step: ' + str(bin_step) + '",\n'
                       '      fontSize: 24\n'
                       '      },\n'
                       '      axisX: {\n'
                       '      title: "Position in genome",\n'
                       '      titleFontSize: 16,\n')
        if indel_flag != []:
            html_out.write('      stripLines:[\n')
            for values in indel_flag[:-1]:
                html_out.write('      {\n'
                               '        startValue: ' + str(values[0]) + ',\n'
                               '        endValue: ' + str(values[1]) + ',\n'
                               '        color: "#16CBEF"\n'
                               '      },\n')
            html_out.write('      {\n'
                           '        startValue: ' + str(indel_flag[-1][0]) + ',\n'
                           '        endValue: ' + str(indel_flag[-1][1]) + ',\n'
                           '        color: "#16CBEF"\n'
                           '      }\n'
                           '      ]\n')
        html_out.write('\n    },\n'
                       '      axisY: {\n'
                       '      title: "Number of reads",\n'
                       '      titleFontSize: 16\n'
                       '},\n'
                       '      data: [')
        for the_data, leg_lab in zip([large_deletions, large_insertions, coverage_array], ['Deletions in read', 'Insertions in read', 'Total reads']):
            html_out.write('         {\n'
                           '         type: "line",\n'
                           '         showInLegend: true,\n'
                           '         legendText: "' + leg_lab + '",\n'
                           '         markerType: "none",\n'
                           '        dataPoints: [\n')
            if refnum == 0:
                bg_out = open(options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() + '.wig', 'w')
                bgt_out = open(options.output_folder + '/bigwig/' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() + '.bwt', 'w')
                bgt_out.write('track type=bigwig bigDataUrl=' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower()
                              + '.bw name=test color=0,0,200 altColor=0,200,0 autoScale=on alwaysZero=on graphType=bar yLineMark=10 yLineOnOff=on\n')
                bgt_out.close()
            else:
                bg_out = open(options.output_folder + '/wiggle/' + leg_lab.replace(' ', '_').replace(')', '').replace('(', '').lower() + '.wig', 'a')
            bg_out.write('fixedStep chrom=' + reference + ' start=1 step=' + str(bin_step) + ' span=' + str(bin_step) + '\n')
            for value in range(0, len(x_axis) - 1):
                html_out.write('{ x: ' + str(x_axis[value]) + ', y: ' + str(the_data[value]) + ' },\n')
                if x_axis[value] >= 0 and x_axis[value] < chrom_size[reference] - bin_step:
                    bg_out.write(str(the_data[value]) + '\n')
            bg_out.close()
            html_out.write('{ x: ' + str(x_axis[-1]) + ', y: ' + str(the_data[-1]) + ' }\n')
            html_out.write('        ]\n      }')
            if not the_data is coverage_array:
                html_out.write(',')
        html_out.write('''      ],
        rangeChanged: syncHandler\n
    });

    chart3.render();
var charts = [chart1, chart2, chart3];

function syncHandler(e) {

    for (var i = 0; i < charts.length; i++) {
        var chart = charts[i];

        if (!chart.options.axisX) chart.options.axisX = {};

        if (e.trigger === "reset") {
            chart.options.axisX.viewportMinimum = chart.options.axisX.viewportMaximum = null;

        } else if (chart !== e.chart) {
            chart.options.axisX.viewportMinimum = e.axisX.viewportMinimum;
            chart.options.axisX.viewportMaximum = e.axisX.viewportMaximum;
        }

        chart.render();

    }
}
function clickHandler(e) {
    var x = parseInt(e.target.id.split(',')[0])
    var y = parseInt(e.target.id.split(',')[1])
    for (var i = 0; i < charts.length; i++) {
        var chart = charts[i];
        chart.options.axisX.viewportMinimum = x;
        chart.options.axisX.viewportMaximum = y;
        chart.render();
        }
}
var zoomButtons = document.getElementsByClassName("zoom");
for (var i = 0; i < zoomButtons.length; i++) {
    var zoomButton = zoomButtons[i]
    zoomButton.addEventListener("click", clickHandler)
    };
  }
  </script>
 <script type="text/javascript" src="/igb/webpage_css_js/canvasjs.min.js"></script></head>
  <h1> ''' + str(html_name) + ''' graphs </h1>
  <div id="chartContainer1" style="height: 400px; width: 100%;">
  </div>
  <br />
    <table style="width: 60%" align="center">
        <tr>
    	    <th>Type</th>
            <th>Position</th>
        </tr>
''')
        for j in clipped_flag:
            html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
            html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
        for j in unclipped_flag:
            html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
            html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
        if unclipped_flag == [] and clipped_flag == []:
            html_out.write('        <tr>\n          <td> no flags </td>\n')
            html_out.write('          <td> no flags </td>\n        </tr>\n')
        html_out.write('''
	</table>
  <br />
  <div id="chartContainer2" style="height: 400px; width: 100%;">
  </div>
  <br />
  <div id="chartContainer3" style="height: 400px; width: 100%;">
  </div>
  <br />
    <table style="width: 60%" align="center">
        <tr>
    	    <th>Type</th>
            <th>Position</th>
        </tr>
''')
        for j in indel_flag:
            html_out.write('        <tr>\n          <td><a class="zoom" id="' + str(j[0] - 1000) + ',' + str(j[1] + 1000) + '"> ' + j[2] + '</a></td>\n')
            html_out.write('          <td>' + str(j[0]) + '..' + str(j[1]) + '</td>\n        </tr>\n')
        if indel_flag == []:
            html_out.write('        <tr>\n          <td> no flags </td>\n')
            html_out.write('          <td> no flags </td>\n        </tr>\n')
        html_out.write('''
	</table>
  <br />
''')
        html_out.write(footer)
        html_out.close()
        if sam.lengths[refnum] >= 300000:
            indexa = 100000/bin_step
        else:
            indexa = sam.lengths[refnum] / 3 / bin_step
        new_array = coverage_array[indexa:-indexa]
        if len(new_array) >= 1000:
            zesteps = len(new_array) / 1000
        else:
            zesteps = 1
        covlist = []
        for jump in range(0, len(new_array) + 1, zesteps):
            covlist.append(sum(new_array[jump:jump+zesteps]) * 1.0/ zesteps)
        out_cov[html_name] = covlist
        out_flag[html_name] = (len(clipped_flag), len(unclipped_flag), len(indel_flag))
    csfile = open(options.output_folder + '/chrom.size', 'w')
    for i in chrom_size:
        csfile.write(i + '\t' + str(chrom_size[i]) + '\n')
    csfile.close()
    for i in os.listdir(options.output_folder + '/wiggle/'):
        subprocess.Popen('wigToBigWig ' + options.output_folder + '/wiggle/' + i + ' ' + options.output_folder + '/chrom.size ' + options.output_folder + '/bigwig/' + i[:-3] + 'bw', shell=True).wait()
    return out_cov, out_flag


# runs BWA on the reference - extends the reference by <max_read_length> on either side to check if contig has been circularised properly
def runBWA(reference, reads, out_dir, num_threads='1', max_read_length=100000):
    if os.path.exists(out_dir) and not os.path.isdir(out_dir):
        sys.stderr.write('Out folder is already a file, please choose another path.')
    elif not os.path.exists(out_dir):
        sys.stdout.write('Creating folder ' + out_dir + '\n')
        os.mkdir(out_dir)
    seqlist = []
    with open(reference) as ref:
        seq = None
        for line in ref:
            if line.startswith('>'):
                if not seq is None:
                    seqlist.append([name, seq])
                name = line.rstrip()
                seq = ''
            else:
                seq += line.rstrip()
    seqlist.append([name, seq])
    out = open(out_dir + '/extended_references.fa', 'w')
    for i in seqlist:
        out.write(i[0] + '\n')
        new_seq = i[1][-max_read_length:] + i[1] + i[1][:max_read_length]
        for j in range(0, len(new_seq), 80):
            out.write(new_seq[j:j+80] + '\n')
    out.close()
    subprocess.Popen('bwa index ' + out_dir + '/extended_references.fa', shell=True).wait()
    subprocess.Popen('bwa mem -t ' + num_threads + ' -k 16 -a ' + out_dir + '/extended_references.fa ' + reads + ' > ' + out_dir + '/alignment.sam', shell=True).wait()
    subprocess.Popen('bwa index ' + reference, shell=True).wait()
    subprocess.Popen('bwa mem -t ' + num_threads + ' -k 16 ' + reference + ' ' + reads + ' > ' + out_dir + '/alignment2.sam', shell=True).wait()
    subprocess.Popen('samtools view -b ' + out_dir + '/alignment2.sam > ' + out_dir + '/alignment.bam', shell=True).wait()
    subprocess.Popen('samtools sort ' + out_dir + '/alignment.bam ' + out_dir + '/alignment.sorted', shell=True).wait()
    subprocess.Popen('samtools index ' + out_dir + '/alignment.sorted.bam', shell=True).wait()


# Create an SVG of blast htis between contigs - also plots coverage on the blast hits
def do_blast(options, header, footer, coverage):
    reference = options.assembly_FASTA
    out_dir = options.output_folder
    subprocess.Popen('makeblastdb -in ' + reference + ' -dbtype nucl -out ' + out_dir + '/all_vs_all_db', shell=True).wait()
    subprocess.Popen('blastn -outfmt 6 -num_threads 6 -query ' + reference + ' -db ' + out_dir + '/all_vs_all_db -out ' + out_dir + '/all_vs_all_blast.out', shell=True).wait()
    blast_dict = {}
    header_names = []
    length_dict = {}
    with open(reference) as ref_file:
        seq = None
        for line in ref_file:
            if line.startswith('>'):
                if not seq is None:
                    header_names.append(name)
                    length_dict[name] = len(seq)
                name = line.rstrip()[1:]
                seq = ''
            else:
                seq += line.rstrip()
        header_names.append(name)
        length_dict[name] = len(seq)
    for i in header_names:
        blast_dict[i] = {}
        for j in header_names:
            blast_dict[i][j] = []
    with open(out_dir + '/all_vs_all_blast.out') as blast_file:
        for line in blast_file:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length = map(int, [qstart, qstop, rstart, rstop, length])
            eval, bitscore, ident = map(float, [eval, bitscore, ident])
            if (length >= 0.1 * length_dict[query] or length >= 10000) and (length != length_dict[subject] or query != subject):
                blast_dict[query][subject].append((qstart, qstop, rstart, rstop, ident))
    width = 3000
    extra = 300
    margin = 50
    y_offset = 125
    hit_thick = 20
    hit_space = 10
    genome_line_thickness = 40
    hit_thick_vert = 5
    svgs_made = []
    ass_name = options.assembly_name
    for i in blast_dict:
        for j in blast_dict[i]:
            if len(blast_dict[i][j]) > 0:
                out_lines1 = []
                out_lines2 = []
                out_lines3 = []
                takenq = [set()]
                takenr = [set()]
                sameq = {}
                samer = {}
                for k in blast_dict[i][j]:
                    if k[2] <= k[3]:
                        fill = (9, 96, 125)
                    else:
                        fill = (125, 22, 36)
                    qx1, qx2 = int(k[0] * float(width)/length_dict[i] + margin), int(k[1] * float(width)/length_dict[i] + margin)
                    if (qx1, qx2) in sameq:
                        row = sameq[(qx1, qx2)]
                        qy = y_offset + row * 20 + 25
                    else:
                        row = None
                        for q, r in enumerate(takenq):
                            hit = False
                            for s in range(qx1 - 10, qx2 + 10):
                                if s in r:
                                    hit = True
                            if not hit:
                                row = q
                                break
                        if row is None:
                            row = len(takenq)
                            takenq.append(set())
                        for s in range(qx1, qx2 + 1):
                            takenq[row].add(s)
                        sameq[(qx1, qx2, fill)] = row
                        qy = y_offset + row * (hit_thick + hit_space) + hit_thick + hit_space + genome_line_thickness
                        out_lines1.append((qx1, qy, qx2, qy, hit_thick, fill))
                    rx1, rx2 = int(k[2] * float(width)/length_dict[j] + margin), int(k[3] * float(width)/length_dict[j] + margin)
                    if rx1 > rx2:
                        rx1, rx2 = rx2, rx1
                    if (rx1, rx2) in samer:
                        row = samer[(rx1, rx2)]
                        ry = y_offset + row * 20 + 25
                    else:
                        row = None
                        for q, r in enumerate(takenr):
                            hit = False
                            for s in range(rx1 - 10, rx2 + 10):
                                if s in r:
                                    hit = True
                            if not hit:
                                row = q
                                break
                        if row is None:
                            row = len(takenr)
                            takenr.append(set())
                        for s in range(rx1, rx2 + 1):
                            takenr[row].add(s)
                        samer[(rx1, rx2, fill)] = row
                        ry = y_offset + row * (hit_thick + hit_space) + hit_space + genome_line_thickness
                        out_lines2.append((rx1, ry, rx2, ry, hit_thick, fill))
                    out_lines3.append(((rx1 + rx2) / 2, ry + hit_thick /2, (qx1 + qx2) / 2, qy + hit_thick / 2, hit_thick_vert, fill))
                height = (len(takenr) + len(takenq)) * (hit_thick + hit_space) + 2 * (y_offset + hit_space + genome_line_thickness) + 100
                svg = scalableVectorGraphics(height, width + 2 * margin + extra)
                for k in out_lines1:
                    svg.drawLine(k[0], k[1], k[2], k[3], k[4], k[5], 0.7)
                for k in out_lines2:
                    svg.drawLine(k[0], height-k[1], k[2], height-k[3], k[4], k[5], 0.7)
                for k in out_lines3:
                    svg.drawLine(k[0], height-k[1], k[2], k[3], k[4], k[5], 0.7)
                x_interval = float(width) / (len(coverage[i[:6]]) -1)
                x_coords = [margin]
                y_coords = [y_offset-genome_line_thickness/2]
                cov_height = 80.0
                cov_max = 1
                for l in coverage[i[:6]]:
                    if l >= cov_max:
                        cov_max = l
                for k,l in enumerate(coverage[i[:6]]):
                    x_coords.append(int(margin + k * x_interval))
                    y_coords.append(y_offset-genome_line_thickness/2-l*1.0/cov_max*cov_height)
                x_coords.append(margin + width)
                y_coords.append(y_offset-genome_line_thickness/2)
                svg.drawPolygon(x_coords, y_coords)
                svg.writeString('Max. coverage', width+2*margin, height/6, 30)
                svg.writeString(str(cov_max), width+2*margin, height/4, 30)
                x_interval = float(width) / (len(coverage[j[:6]]) -1)
                x_coords = [margin]
                y_coords = [height-y_offset+genome_line_thickness/2]
                cov_height = 80.0
                cov_max = 1
                for l in coverage[j[:6]]:
                    if l >= cov_max:
                        cov_max = l
                for k,l in enumerate(coverage[j[:6]]):
                    x_coords.append(int(margin + k * x_interval))
                    y_coords.append(height - y_offset + genome_line_thickness/2+l*1.0/cov_max*cov_height)
                x_coords.append(margin + width)
                y_coords.append(height- y_offset+genome_line_thickness/2)
                svg.drawPolygon(x_coords, y_coords)
                svg.writeString('Max. coverage', width+2*margin, height*5/6, 30)
                svg.writeString(str(cov_max), width+2*margin, height*11/12, 30)
                svg.drawLine(margin,y_offset, width + margin, y_offset, genome_line_thickness)
                svg.writeString(i + ' length=' + str(length_dict[i]), margin, y_offset + genome_line_thickness/4, 30)
                svg.drawLine(margin,height-y_offset, width+ margin, height-y_offset, genome_line_thickness)
                svg.writeString(j + ' length=' + str(length_dict[j]), margin, height-y_offset + genome_line_thickness/4, 30)
                svg.drawLine(width+2*margin,height*3/8, width+2*margin+40, height*3/8, hit_thick, (9, 96, 125))
                svg.writeString('Direct hit', width+2*margin + 50, height*3/8+ genome_line_thickness/4, 30)
                svg.drawLine(width+2*margin, height*5/8, width+2*margin+40, height*5/8, hit_thick, (125, 22, 36))
                svg.writeString('Inverted hit', width+2*margin + 50, height*5/8+ genome_line_thickness/4, 30)
                svg.writesvg(out_dir + '/qc_website/blast/' + i + '_' + j + '.svg')
                svgs_made.append('/igb/' + ass_name +  '/qc_website/blast/' + i + '_' + j + '.svg')
    for i in header_names:
        html_name = i[:6]
        html_out = open(out_dir + '/qc_website/blast/' + html_name + '_blast.html', 'w')
        html_out.write(header)
        html_out.write('<h1> ' + html_name + ' blast results </h1>\n')
        for j in svgs_made:
            if i in j:
                html_out.write('<img src="' + j + '" alt="" width="100%" />')
        html_out.write(footer)



# creates the header/footer sidebar for each page so they have a consistent theme
def get_page_bookends(options):
    reference = options.assembly_FASTA
    ass_name = options.assembly_name
    namelist = []
    with open(reference) as ref:
        for line in ref:
            if line.startswith('>'):
                namelist.append(line.rstrip()[1:])
    namelist.sort()
    header = '''<!--Force IE6 into quirks mode with this comment tag-->
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<title>Pacbio assembly QC</title>
<link rel="stylesheet" type="text/css" href="/igb/webpage_css_js/template.css" />
<script src="/igb/webpage_css_js/svg-pan-zoom.js"></script>
</head>

<body>

<div id="framecontent">
<div class="innertube">
'''
    try:
        header += '<h1>Assembly ' + ass_name.split('_')[2] + ' TOC</h1>\n'
    except IndexError:
        header += '<h1>Assembly ' + ass_name + ' TOC</h1>\n'
    header += '<h3><a href="/igb/' + ass_name + '/index.html">Overview</a></h3>\n'
    for i in namelist:
        num = i[:6]
        header += '<h3> ' + num + '</h3>\n'
        header += '<p><a href="/igb/' + ass_name + '/qc_website/graphs/' + num + '_graphs.html">Graphs</a></p>\n'
        header += '<p><a href="/igb/' + ass_name + '/qc_website/blast/' + num + '_blast.html">BLAST</a></p>\n'
    header += '<h3>Downloads</h3>\n'
    header += '<p><a href="/igb/' + ass_name + '/' + ass_name + '.fasta"> FASTA </ a></p>\n'
    header += '<p><a href="/igb/' + ass_name + '/' + ass_name + '.bed"> bed </a></p>\n'
    header += '<p><a href="/igb/' + ass_name + '/' + ass_name + '.2bit"> 2bit </a></p>\n'
    header += '<p><a href="/igb/' + ass_name + '/annots.xml"> Annotation </a></p>\n'
    header += '<p><a href="/igb/' + ass_name + '/genome.txt"> Genome </a></p>\n'
    header += '<p><a href="/igb/' + ass_name + '/mlst.txt"> MLST </a></p>\n'
    header += '''

</div>
</div>

<div id="maincontent">
<div class="innertube">

'''
    footer = '''
<div style="margin-top: 1em">
</div>

<div style="text-align: center; font: normal 6px Arial; margin-top: 1em;"> This page was generated by createQCwebpage.py - for help email: Mitchell Sullivan (mjsull@gmail.com) </a></div>

</div>
</div>


</body>
</html>'''
    return header, footer


# Creates a SVG of the assembly graph file
def create_graph(options):
    gkp_store = options.gkp_location
    tig_store = options.tig_location
    best_edge = options.best_edges
    output = options.output_folder + '/qc_website/graph.svg'
    G=nx.DiGraph()
    frg_to_tig = {}
    cout = {}
    args = shlex.split("tigStore -g %s -t %s 2 -D unitiglist" % (gkp_store, tig_store ))
    out = subprocess.check_output(args)
    out = out.split("\n")
    for l in out:
        l = l.strip().split()
        if len(l) == 0: continue
        if l[0] == "maID": continue
        unitig_id = int(l[0])
        os.system("tigStore -g %s -t %s 2 -d frags -u %d > frag_list" % ( gkp_store, tig_store, unitig_id) )
        args = shlex.split( "tigStore -g %s -t %s 2 -d frags -u %d" % ( gkp_store, tig_store, unitig_id) )
        f_out = subprocess.check_output(args)
        f_out = f_out.split("\n")
        for l in f_out:
            """FRG    1453 179419,182165"""
            l = l.replace(",", " ")
            l = l.strip().split()
            if len(l) == 0: continue
            frg_id = l[1]
            frg_to_tig[frg_id] = unitig_id
    with open(best_edge) as f:
        for l in f:
            if l[0] == "#": continue
            l = l.strip().split()
            id1, lib_id, best5, o1, best3, o3, j1, j2 = l
    #        id1, lib_id, best5, o1, best3, o3 = l
            try:
                G.add_node(id1, label="utg%s" % frg_to_tig[id1], size=int(cout["unitig_%s"%frg_to_tig[id1]][0]), cov=float(cout["unitig_%s"%frg_to_tig[id1]][1]),ref=(cout["unitig_%s"%frg_to_tig[id1]][2]))
            except KeyError:
                G.add_node(id1, label="utg%s" % frg_to_tig[id1], size=int(0), cov=float(0))
            if best5 != "0":
                G.add_edge(best5, id1)
            if best3 != "0":
                G.add_edge(id1, best3)
    Agraph = nx.nx_agraph.to_agraph(G)
    Agraph.edge_attr['penwidth'] = 5
    Agraph.node_attr['style'] = 'filled'
    h, s, l = 175, 0.8, 0.5
    color_dict = {}
    for nodes in Agraph.nodes_iter():
        n = Agraph.get_node(nodes)
        label = n.attr['label']
        if not label in color_dict:
            h = (h + 65) % 360
            color_dict[label] = hsl_to_colorstr(h, s, l)
        color = color_dict[label]
        n.attr['fillcolor']=color
    Agraph.draw(output, prog='neato', format='svg')

# creates a dot plot of the assembly
def create_dot_plot(options):
    reference = options.assembly_FASTA
    header_names = []
    min_length = 500
    min_ident = 95
    out_dir = options.output_folder
    with open(reference) as ref_file:
        seq = None
        for line in ref_file:
            if line.startswith('>'):
                if not seq is None:
                    header_names.append((len(seq), name))
                name = line.rstrip()[1:]
                seq = ''
            else:
                seq += line.rstrip()
        header_names.append((len(seq), name))
    draw_list = []
    with open(out_dir + '/all_vs_all_blast.out') as blast_file:
        for line in blast_file:
            query, subject, ident, length, mm, indel, qstart, qstop, rstart, rstop, eval, bitscore = line.split()
            qstart, qstop, rstart, rstop, length = map(int, [qstart, qstop, rstart, rstop, length])
            eval, bitscore, ident = map(float, [eval, bitscore, ident])
            if length >= min_length and ident >= min_ident:
                draw_list.append((query, subject, qstart, qstop, rstart, rstop))
    gap_bp = 50000
    size = 5000
    pad = 100
    header_names.sort(reverse=True)
    curr_x = 0
    pos_dict = {}
    for i in header_names:
        pos_dict[i[1]] = curr_x
        curr_x += i[0] + gap_bp

    total_length = float(curr_x - gap_bp)
    svg = scalableVectorGraphics(size + 2*pad + 1000, size + 2*pad+ 2000)
    curr_x = 0
    for i in header_names:
        curr_y = total_length
        for j in header_names:
            curr_y -= j[0]
            svg.drawOutRect(curr_x / total_length * size + pad, curr_y / total_length * size + pad, i[0] / total_length * size, j[0] / total_length * size, (255, 255, 255), 5, 0)
            curr_y -= gap_bp
        svg.writeString(i[1], (curr_x + i[0]/2) /total_length * size + pad-24, pad + pad/2 + size, 96, rotate=-1)
        svg.writeString(i[1], pad + pad/2 + size, size - ((curr_x + i[0]/2) /total_length * size) + pad+24, 96)
        curr_x += i[0] + gap_bp
    for i in draw_list:
        query, subject, qstart, qstop, rstart, rstop = i
        svg.drawLine((pos_dict[subject] + rstart) / total_length * size + pad,
                     (size - (pos_dict[query] + qstart) / total_length * size) + pad,
                     (pos_dict[subject] + rstop) / total_length * size + pad,
                     (size - (pos_dict[query] + qstop) / total_length * size) + pad, 5)
    svg.writesvg(out_dir + '/qc_website/dot_plot.svg')



# create the index page - contains table of information, dotplot and graph
def write_index(options, header, footer, coverage, flags):
    contigs = []
    contig_names = []
    contig_lengths = {}
    with open(options.assembly_FASTA) as assem_file:
        seq = None
        for line in assem_file:
            if line.startswith('>'):
                if not seq is None:
                    contig_names.append(name)
                    contig_lengths[name] = str(len(seq))
                name = line.rstrip()[1:]
                seq = ''
            else:
                seq += line.rstrip()
    contig_names.append(name)
    contig_lengths[name] = str(len(seq))
    assign_dict = {}
    if os.path.exists(options.output_folder + '/assignments.tsv'):
        with open(options.output_folder + '/assignments.tsv') as assignments:
            for line in assignments:
                assign_dict[line.split('\t')[0]] = line.split('\t')[1]
    else:
        for q in contig_names:
            assign_dict[q] = 'none'
    for q in contig_names:
        if q[6] == 'x':
            trimmed = 'n'
        else:
            trimmed = 'y'
        if q[7] == 'x':
            reorientated = 'n'
        else:
            reorientated = 'y'
        if q[8] == 'x':
            ends_corrected = 'n'
        else:
            ends_corrected = 'y'
        contigs.append((q, assign_dict[q], contig_lengths[q], sum(coverage[q[:6]]) * 1.0 / len(coverage[q[:6]]), trimmed, ends_corrected, reorientated, flags[q[:6]][0], flags[q[:6]][1], flags[q[:6]][2]))
    html_out = open(options.output_folder + '/index.html', 'w')
    html_out.write(header)
    html_out.write('''
  <br />
  <h2> Summary of contigs </h2>
    <table style="width: 90%" align="center">
        <tr>
    	    <th>Contig Name</th>
            <th>Designation</th>
            <th>Length</th>
            <th>Mean coverage</th>
            <th>Circularised</th>
            <th>Polished</th>
            <th>Reorientated</th>
            <th>Clipped read flagged</th>
            <th>Unclipped flagged</th>
            <th>Indel flagged</th>
        </tr>
''')
    for i in contigs:
        html_out.write('        <tr>\n')
        for j in i:
            html_out.write('          <td>' + str(j) + '</td>\n')
        html_out.write('        </tr>\n')
    html_out.write('''
	</table>
	<br >
  <h2> Dot Plot </h2>
    <img src="qc_website/dot_plot.svg" alt="" width="100%" />
  <h2>Assembly graph</h2>
    <object id="demo-tiger" type="image/svg+xml" data="qc_website/graph.svg" style="width: 50%; height: 50%; border:1px solid black; ">Your browser does not support SVG</object>
    <script>
      window.onload = function() {
        svgPanZoom('#demo-tiger', {
          zoomEnabled: true,
          controlIconsEnabled: true,
          maxZoom: 100
        });
      };
    </script>


  <br />
  <h2> Log of changes to assembly file. </h2>
''')
    if not os.path.exists(options.output_folder + '/log.tsv'):
        log = open(options.output_folder + '/log.tsv', 'w')
        log.write(str(datetime.datetime.now()).split('.')[0] + '\tQC webpage created\tnone')
        log.close()
    with open(options.output_folder + '/log.tsv') as log:
        html_out.write('\n    <table style="width: 90%" align="center">\n'
                       '        <tr>\n'
                       '    	    <th>Time</th>\n'
                       '            <th>Event</th>\n'
                       '            <th>Notes</th>\n'
                       '        </tr>\n')
        for line in log:
            time, event, note = line.split('\t')
            html_out.write('    <tr>\n'
                           '        <td>' + time + '</td>\n'
                           '        <td>' + event + '</td>\n'
                           '        <td>' + note + '</td>\n    </tr>\n')
            html_out.write('<!-- Table ends here -->\n')
            html_out.write('</table>\n')


    html_out.write(footer)


# Initially this was going to be a way to log manual changes to the assembly - may implement at some stage..
def log_change(options):
    out_html = ''
    time = str(datetime.datetime.now()).split('.')[0]
    if options.log_change_note is None:
        note = 'none'
    else:
        note = options.log_change_note
    event = options.log_change
    with open(options.output_folder + '/index.html') as index:
        for line in index:
            if line == '<!-- Table ends here -->\n':
                out_html += ('    <tr>\n'
                           '        <td>' + time + '</td>\n'
                           '        <td>' + event + '</td>\n'
                           '        <td>' + note + '</td>\n    </tr>\n<!-- Table ends here -->\n')
            else:
                out_html += line
    out = open(options.output_folder + '/index.html', 'w')
    out.write(out_html)

# copy all the required html files to html_loc
def copy_html(options):
    if not os.path.exists(options.html_loc):
        os.makedirs(options.html_loc)
    shutil.copy(options.output_folder + '/index.html', options.html_loc + '/index.html')
    if os.path.exists(options.html_loc + '/qc_website/'):
        shutil.rmtree(options.html_loc + '/qc_website/')
    try:
        shutil.copytree(options.output_folder + '/qc_website/', options.html_loc + '/qc_website/')
    # Directories are the same
    except shutil.Error as e:
        sys.sterr.write('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        sys.sterr.write('Directory not copied. Error: %s' % e)
# copy the bedgraph files to html_loc
def copy_bedgraph(options):
    if os.path.exists(options.html_loc + '/wiggle/'):
        shutil.rmtree(options.html_loc + '/wiggle/')
    try:
        shutil.copytree(options.output_folder + '/wiggle/', options.html_loc + '/wiggle/')
    # Directories are the same
    except shutil.Error as e:
        sys.sterr.write('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        sys.sterr.write('Directory not copied. Error: %s' % e)


parser = OptionParser()
parser.add_option("-o", "--output_folder", help="write results to FOLDER", metavar="FOLDER")
parser.add_option("-w", "--html_loc", help="Location to copy files for webpage.", metavar="/directory/")
parser.add_option("-b", "--bin_size", type="int", default=1000, help="Bin size for drawing graph")
parser.add_option("-s", "--bin_step", type="int", default=100, help="Size of bin steps for drawing graph.")
parser.add_option("-i", "--indel_min", type="int", default=10, help="minimum indel size to flag in indel graph")
parser.add_option("-x", "--flag_ratio", type="float", default=0.5, help="minimum indel size to flag in indel graph")
parser.add_option("-f", "--assembly_FASTA", help="Run create_QC_webpage on FASTA, will try determine assembly name from filename. Need to provide script with read file as well", metavar="ASSEMBLY_FASTA")
parser.add_option("-r", "--assembly_reads", help="Read file for mapping to ASSEMBLY_FASTA", metavar="Reads.fa")
parser.add_option("-g", "--graph_dir", help="Directory with best.edges, tigStore and gkpStore")
parser.add_option("-a", "--assembly_name", help="Assembly name i.e. C_difficile_CD00001_1A_123456")
(options, args) = parser.parse_args()

options.gkp_location = options.graph_dir + '/celera-assembler.gkpStore'
options.tig_location = options.graph_dir + '/celera-assembler.tigStore'
if os.path.exists(options.graph_dir + '/best.edges'):
    options.best_edges = options.graph_dir + '/best.edges'
else:
    options.best_edges = options.graph_dir + '/4-unitigger/best.edges'


if not os.path.exists(options.assembly_FASTA):
    sys.stderr.write('We expected the assembly FASTA to be here: ' + options.assembly_FASTA + ' but no file was found. \n')
    sys.exit()
if not os.path.exists(options.assembly_reads):
    sys.stderr.write('We expected the read FASTQ to be here: ' + options.assembly_reads + ' but no file was found.\n')
    sys.exit()

large_gap_min = options.indel_min
if not os.path.exists(options.output_folder + '/qc_website'):
    os.makedirs(options.output_folder + '/qc_website')
    os.makedirs(options.output_folder + '/qc_website/blast')
    os.makedirs(options.output_folder + '/qc_website/graphs')
else:
    if os.path.isdir(options.output_folder + '/qc_website'):
        if not os.path.exists(options.output_folder + '/qc_website/blast'):
            os.makedirs(options.output_folder + '/qc_website/blast')
        if not os.path.exists(options.output_folder + '/qc_website/graphs'):
            os.makedirs(options.output_folder + '/qc_website/graphs')
    else:
        sys.stderr.write('File at qc_website, folder could not be created. Exiting...')
        sys.exit()
if not os.path.exists(options.output_folder + '/wiggle'):
    os.makedirs(options.output_folder + '/wiggle')
elif not os.path.isdir(options.output_folder + '/wiggle'):
    sys.stderr.write('File at wiggle, folder could not be created. Exiting...')
    sys.exit()
if not os.path.exists(options.output_folder + '/bigwig'):
    os.makedirs(options.output_folder + '/bigwig')
elif not os.path.isdir(options.output_folder + '/bigwig'):
    sys.stderr.write('File at bigwig, folder could not be created. Exiting...')
    sys.exit()

header, footer = get_page_bookends(options)
runBWA(options.assembly_FASTA, options.assembly_reads, options.output_folder, '3')
coverage, out_flag = draw_graph(options, header, footer)
do_blast(options, header, footer, coverage)
create_dot_plot(options)
create_graph(options)
write_index(options, header, footer, coverage, out_flag)
if not options.html_loc is None:
    copy_html(options)
    #copy_bedgraph(options)
