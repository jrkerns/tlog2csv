__author__ = 'James Kerns'
__version__ = '0.1.0'

"""An executable script that converts a single Varian TrueBeam trajectory log
(*.bin file) to a comma-separated variable (CSV) file, similar to dynalog CSV files.

Algorithm implementation is based on the log reader of Pylinac.


Copyright (c) 2014 James Kerns

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without
limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions
of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

"""

import csv
import struct
import os.path as osp
import sys
try:
    # Python 3.x
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename, askopenfilenames, askdirectory
except ImportError:
    # Python 2.x
    from Tkinter import Tk
    from tkFileDialog import askopenfilename, askopenfilenames, askdirectory


class Axis(object):
    def __init__(self, actual, expected=None):
        self.actual = actual
        if expected is not None:
            self.expected = expected

def decode_binary(filecontents, dtype, num_values=1, cursor_shift=0):
    """This method is the main "decoder" for reading in trajectory log binary data into human data types.

    :param filecontents: the complete file having been read with .read().
    :param dtype: the expected data type to return. If int or float, will return numpy array.
    :type dtype: str, int, float
    :param num_values: the expected number of dtype to return; note that this is not the same as the # of bytes.
    :type num_values: int
    :param cursor_shift: the number of bytes to move the cursor forward after decoding. This is used if there is a
        reserved section after the read-in segment.
    :type cursor_shift: int
    """
    global cursor
    fc = filecontents

    if dtype == str:  # if string
        output = fc[cursor:cursor + num_values]
        if type(fc) is not str:  # in py3 fc will be bytes
            output = output.decode()
        # Strip the padding ("\x00")
        output = output.strip('\x00')
        cursor += num_values
    elif dtype == int:
        ssize = struct.calcsize('i') * num_values
        output = struct.unpack('i' * num_values, fc[cursor:cursor + ssize])
        if len(output) == 1:
            output = int(output[0])
        cursor += ssize
    elif dtype == float:
        ssize = struct.calcsize('f') * num_values
        output = struct.unpack('f' * num_values, fc[cursor:cursor + ssize])
        if len(output) == 1:
            output = float(output[0])
        cursor += ssize
    else:
        raise TypeError("decode_binary datatype was not valid")

    cursor += cursor_shift  # shift cursor if need be (e.g. if a reserved section follows)
    return output


"""Get the Tlog file path"""
# if executing stand-alone, use a UI dialog box to get the
# Tlog file.
if len(sys.argv) == 1:
    Tk().withdraw()
    tlog_file = askopenfilename()
else:
    raise NotImplementedError("Command line use not implemented yet. Please use standalone.")
    #TODO: add command line arguments
if not osp.isfile(tlog_file):
    raise IOError("Target specified was not a file")

"""Read in Tlog"""
# read in trajectory log binary data
fcontent = open(tlog_file, 'rb').read()

cursor = 0

# Unpack the content according to respective section and data type (see log specification file).
try:
    signature = decode_binary(fcontent, str, 16)  # for version 1.5 will be "VOSTL"
except:
    raise IOError("Target specified does not seem to be a trajectory log file.")
version = float(decode_binary(fcontent, str, 16))  # in the format of 2.x or 3.x
header_size = decode_binary(fcontent, int)  # fixed at 1024 in 1.5 specs
sampling_interval = decode_binary(fcontent, int)
num_axes = decode_binary(fcontent, int)
axis_enum = decode_binary(fcontent, int, num_axes)
samples_per_axis = decode_binary(fcontent, int, num_axes)
num_mlc_leaves = samples_per_axis[-1] - 2  # subtract 2 (each carriage counts as an "axis" and must be removed)
# cursor = num_axes * 4  # there is a reserved section after samples per axis. this moves it past it.
clinac_scale = decode_binary(fcontent, int)
num_subbeams = decode_binary(fcontent, int)
is_truncated = decode_binary(fcontent, int)
num_snapshots = decode_binary(fcontent, int)
# the section after MLC model is reserved. Cursor is moved to the end of this reserved section.
mlc_model = decode_binary(fcontent, int, cursor_shift=1024 - (64 + num_axes * 8))

# read in subbeam data. These are for auto-sequenced beams. If not autosequenced, separate logs are created.
# Currently there is no good way of dealing with this data, but fortunately autosequencing is rare at this time.
if num_subbeams:
    for beam in range(num_subbeams):
        cursor += 80

# ----------------------------------------------------------------------
# assignment of snapshot data (actual & expected of MLC, Jaw, Coll, etc)
#----------------------------------------------------------------------

# step size in bytes
step_size = sum(samples_per_axis) * 2

# read in all snapshot data at once, then assign
snapshot_data = decode_binary(fcontent, float, step_size * num_snapshots)

# collimator
collimator = Axis(expected=snapshot_data[0::step_size],
                       actual=snapshot_data[1::step_size])

# gantry
gantry = Axis(expected=snapshot_data[2::step_size],
                   actual=snapshot_data[3::step_size])

# jaws
jaw_y1 = Axis(expected=snapshot_data[4::step_size],
                   actual=snapshot_data[5::step_size])
jaw_y2 = Axis(expected=snapshot_data[6::step_size],
                   actual=snapshot_data[7::step_size])
jaw_x1 = Axis(expected=snapshot_data[8::step_size],
                   actual=snapshot_data[9::step_size])
jaw_x2 = Axis(expected=snapshot_data[10::step_size],
                   actual=snapshot_data[11::step_size])

# couch
couch_vrt = Axis(expected=snapshot_data[12::step_size],
                      actual=snapshot_data[13::step_size])
couch_lng = Axis(expected=snapshot_data[14::step_size],
                      actual=snapshot_data[15::step_size])
couch_lat = Axis(expected=snapshot_data[16::step_size],
                      actual=snapshot_data[17::step_size])
couch_rtn = Axis(expected=snapshot_data[18::step_size],
                      actual=snapshot_data[19::step_size])

# MU
mu = Axis(expected=snapshot_data[20::step_size],
               actual=snapshot_data[21::step_size])

# beam hold state
beam_hold = Axis(expected=snapshot_data[22::step_size],
                      actual=snapshot_data[23::step_size])

# control point
control_point = Axis(expected=snapshot_data[24::step_size],
                          actual=snapshot_data[25::step_size])

# carriages
carraigeA = Axis(expected=snapshot_data[26::step_size],
                      actual=snapshot_data[27::step_size])
carriageB = Axis(expected=snapshot_data[28::step_size],
                      actual=snapshot_data[29::step_size])

# assign MLC data for all leaves. Usually 1-60 is bank A, 61-120 is bank B
# Units are in cm.
mlc = {}
for leaf in range(num_mlc_leaves):
    mlc[leaf] = (Axis(expected=snapshot_data[30 + (2 * leaf)::step_size],
                          actual=snapshot_data[31 + (2 * leaf)::step_size]))


def write_single_value(writer, description, value, unit=None):
    writer.writerow([description, str(value), unit])

def write_array(writer, description, value, unit=None):
        # write expected
        if unit is None:
            writer.writerow([description + ' Expected'])
        else:
            writer.writerow([description + ' Expected in units of', unit])
        writer.writerow(value.expected)
        # write actual
        if unit is None:
            writer.writerow([description + ' Actual'])
        else:
            writer.writerow([description + ' Actual in units of ', unit])
        writer.writerow(value.actual)


"""Write the data to CSV"""
csv_filename = tlog_file.replace('bin', 'csv')
with open(csv_filename, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    write_single_value(writer, 'Tlog File:', tlog_file)
    write_single_value(writer, 'Signature:', signature)
    write_single_value(writer, 'Version:', version)
    write_single_value(writer, 'Header Size:', header_size)
    write_single_value(writer, 'Sampling Inteval:', sampling_interval, 'ms')
    write_single_value(writer, 'Number of Axes:', num_axes)
    write_single_value(writer, 'Axis Enumeration:', axis_enum)
    write_single_value(writer, 'Samples per Axis:', samples_per_axis)
    write_single_value(writer, 'Axis Scale:', clinac_scale)
    write_single_value(writer, 'Number of Subbeams:', num_subbeams)
    write_single_value(writer, 'Is Truncated?', is_truncated)
    write_single_value(writer, 'Number of Snapshots:', num_snapshots)
    write_single_value(writer, 'MLC Model:', mlc_model)
    write_array(writer, 'Gantry', gantry, 'degrees')
    write_array(writer, 'Collimator', collimator, 'degrees')
    write_array(writer, 'Couch Lat', couch_lat, 'cm')
    write_array(writer, 'Couch Lng', couch_lng, 'cm')
    write_array(writer, 'Couch Rtn', couch_rtn, 'degrees')
    write_array(writer, 'MU', mu, 'MU')
    write_array(writer, 'Beam Hold', beam_hold)
    write_array(writer, 'Control Point', control_point)
    write_array(writer, 'Carriage A', carraigeA, 'cm')
    write_array(writer, 'Carriage B', carriageB, 'cm')
    for leaf in mlc:
        write_array(writer, 'Leaf ' + str(leaf+1), mlc[leaf], 'cm')


