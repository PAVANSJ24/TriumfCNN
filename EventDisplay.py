#!/usr/bin/env python
# coding: utf-8
# # EventDisplay -- to display mPMT events in new npz file format
# Edit to input the full geometry file, and npz data file that your are interested in.
# Authors: Blair Jamieson, Connor Boubard
# June 2020

# In[138]:
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math
import random

# In[2]:
datafile = np.load('event998.npz', allow_pickle=True)
geofile = np.load('mpmt_full_geo.npz', allow_pickle=True)
# # First let's explore the geometry file
# Make sure we can find the phototube locations, and build a mapping from the three dimensional locations of the PMTs.


# In[3]:
geofile.files

# In[4]:
tubes = geofile['tube_no']
tubes

# In[5]:
tube_xyz = geofile['position']
tube_x = tube_xyz[:, 0]
tube_y = tube_xyz[:, 1]
tube_z = tube_xyz[:, 2]
R = (tube_x.max() - tube_x.min()) / 2.0
H = tube_y.max() - tube_y.min()
print("R=", R, "H=", H)
print("min_x=", tube_x.min(), "max_x=", tube_x.max(), "diameter=", tube_x.max() - tube_x.min())
print("min_z=", tube_z.min(), "max_z=", tube_z.max(), "diameter=", tube_z.max() - tube_z.min())
print("min_y=", tube_y.min(), "max_y=", tube_y.max(), "height=", tube_y.max() - tube_y.min())

# In[6]:
tube_dir = geofile['orientation']

# In[7]:
fig = plt.figure(figsize=[15, 15])
ax = fig.add_subplot(111, projection='3d')
ax.scatter(tube_x, tube_y, tube_z, marker='.')
ax.set_xlabel('X (cm)')
ax.set_ylabel('Y (cm)')
ax.set_zlabel('Z (cm)')
ax.view_init(elev=45.0, azim=45.0)
plt.show()


# In[8]:
def PMT_to_flat_cylinder_mapping(tubes, tube_xyz):
    """
    Build dictionary of PMT number, to (x,y) on a flat cylinder
    
    N.B. Tube numbers in full geometry file go from 1:NPMTs, but it seems like
    the event data number from 0:NPMTs-1, so subtracting 1 from tube number here?
    """
    mapping = {}
    for idx, tube in enumerate(tubes):
        x = tube_xyz[idx, 0]
        y = tube_xyz[idx, 1]
        z = tube_xyz[idx, 2]
        if y > 500.:
            # in top circle of cylinder
            xflat = x
            yflat = 897.6 + z
            mapping[int(tube - 1)] = [xflat, yflat]

        elif y < -500.:
            # in bottom circle of cylinder
            xflat = x
            yflat = -897.6 + z
            mapping[int(tube - 1)] = [xflat, yflat]

        else:
            # in barrel part of cylinder
            theta = math.atan2(z, x)
            xflat = R * theta
            yflat = y
            mapping[int(tube - 1)] = [xflat, yflat]
    return mapping


PMTFlatMapping = PMT_to_flat_cylinder_mapping(tubes, tube_xyz)


# In[9]:


def PMT_to_flat_cylinder_map_positive(tubes, tube_xyz):
    """
    Build dictionary of PMT number, to (x,y) on a flat cylinder
    
    N.B. Tube numbers in full geometry file go from 1:NPMTs, but it seems like
    the event data number from 0:NPMTs-1, so subtracting 1 from tube number here?
    
    """
    mapping = {}
    for idx, tube in enumerate(tubes):
        x = tube_xyz[idx, 0]
        y = tube_xyz[idx, 1]
        z = tube_xyz[idx, 2]
        if y > 500.:
            # in top circle of cylinder
            xflat = x + 1162.7
            yflat = 2165.2 + z
            mapping[int(tube - 1)] = [int(round(xflat)), int(round(yflat))]

        elif y < -500.:
            # in bottom circle of cylinder
            xflat = x + 1162.7
            yflat = 370.1 + z
            mapping[int(tube - 1)] = [int(round(xflat)), int(round(yflat))]

        else:
            # in barrel part of cylinder
            theta = math.atan2(z, x)
            xflat = R * theta + 1162.7
            yflat = y + 1267.7
            mapping[int(tube - 1)] = [int(round(xflat)), int(round(yflat))]
    return mapping


PMTFlatMapPositive = PMT_to_flat_cylinder_map_positive(tubes, tube_xyz)
# PMTFlatMapPositive


# In[10]:


xflatvals = []
yflatvals = []
for tube in PMTFlatMapping:
    xflatvals.append(PMTFlatMapping[tube][0])
    yflatvals.append(PMTFlatMapping[tube][1])

fig = plt.figure(figsize=[50, 50])
plt.plot(xflatvals, yflatvals, '.', color='black')
# uncomment the following line to save the image to a pdf file
# plt.savefig( 'mpmt_flatmap_positions.pdf')


# Try making plot of PMT locations as image map

# In[11]:


fig = plt.figure(figsize=[40, 40])

preimage = np.zeros([2506, 2317])
for tube in PMTFlatMapPositive:
    for dx in range(-3, 4):
        for dy in range(-3, 4):
            if abs(dx) == 3 and abs(dy) == 3:
                continue
            preimage[PMTFlatMapPositive[tube][1] + dx, PMTFlatMapPositive[tube][0] + dy] = tube + 5000
plt.imshow(preimage)
fig.suptitle('PMT Tube Number + 5000', fontsize=30)
plt.xlabel('Distance CCW on perimeter from x-axis (cm)', fontsize=24)
plt.ylabel('Y (cm)', fontsize=24)
# plt.set_cmap('Greys')
plt.set_cmap('hot_r')
# plt.set_cmap('gist_heat_r')
plt.colorbar()
# plt.savefig( 'mpmt_flatmap_positions_img.pdf',dpi=300)


# # Now lets look at a few events
# 
# Need to explore contents of event data, and make presentations of that data.  Things to present are:
# 
# - Time of hits on each PMT
# - Charge of hits on each PMT

# In[12]:


datafile.files

# In[40]:


evno = 0


def GetEvent(event_no):
    dtubes = datafile['digi_hit_pmt'][event_no]
    dcharges = datafile['digi_hit_charge'][event_no]
    dtimes = datafile['digi_hit_time'][event_no]
    return (dtubes, dcharges, dtimes)


digitubes, digicharges, digitimes = GetEvent(evno)

# In[14]:


number_of_events = len(datafile['digi_hit_pmt'])
number_of_events

# In[15]:


digicharges

# In[16]:


digitubes

# In[17]:


len(digitubes)

# In[18]:


len(digicharges)


# In[19]:


def EventDisplay(tubes, quantities, title="Charge", cutrange=[-1, -1]):
    """
    tubes == np.array of PMTs that were hit
    quantities == np.array of PMT quantities (either charge or time)
    title == title to add to display
    cutrange == minimum and maximum values on plot (or set both same for default)
    """

    fig = plt.figure(figsize=[12, 12])
    preimage = np.zeros([2506, 2317])
    # maxquantity = quantities.max()
    # preimage *= maxquantity*1.2
    imgmin = quantities.min()
    imgmax = quantities.max()
    for idx, tube in enumerate(tubes):
        if cutrange[0] != cutrange[1]:
            if quantities[idx] < cutrange[0] or quantities[idx] > cutrange[1]:
                continue
        for dx in range(-3, 4):
            for dy in range(-3, 4):
                if abs(dx) == 3 and abs(dy) == 3:
                    continue

                # print( "idx=", idx, " len(quantities)=",len(quantities), " tube=", tube, " len(PMTFlatMap)=", len(PMTFlatMapPositive))
                preimage[PMTFlatMapPositive[tube][1] + dx, PMTFlatMapPositive[tube][0] + dy] = quantities[idx]

    if cutrange[0] != cutrange[1]:
        imgmin = cutrange[0]
        imgmax = cutrange[1]
    plt.imshow(preimage, extent=[-1162.7, 1162.7, -1267.7, 1267.7], vmin=imgmin, vmax=imgmax)
    fig.suptitle(title, fontsize=20)
    plt.xlabel('Distance CCW on perimeter from x-axis (cm)', fontsize=18)
    plt.ylabel('Y (cm)', fontsize=16)
    # plt.set_cmap('YlGnBu')
    plt.set_cmap('cubehelix_r')
    # plt.set_cmap('gnuplot2_r')
    # plt.set_cmap('gist_heat_r')
    # plt.set_cmap('inferno_r')
    # plt.set_cmap('pink_r')
    plt.colorbar()


# In[20]:


EventDisplay(digitubes, digicharges)
# uncomment following line to make pdf image
# plt.savefig( 'event0_charge.pdf',dpi=300)


# In[21]:


EventDisplay(digitubes, digitimes, "Time")


# plt.savefig( 'event0_time.pdf',dpi=300)


# # Make histograms of the charge and time distributions
# 
# 

# In[22]:


def EventDisplayHist(quantities, title="Charge", cutrange=[-1, -1]):
    """
    quantities = np.array of values to histogram
    title = title to add to the histogram
    cutrange = x-axis range to include in histogram
    Makes a histogram with 100 bins
    """
    fig = plt.figure(figsize=[12, 12])
    imgmin = quantities.min()
    imgmax = quantities.max()
    if cutrange[0] != cutrange[1]:
        imgmin = cutrange[0]
        imgmax = cutrange[1]
    plt.hist(quantities, 100, [imgmin, imgmax])
    # fig.suptitle(title, fontsize=20)
    plt.xlabel(title, fontsize=18)
    plt.ylabel('Count / bin', fontsize=16)


EventDisplayHist(digicharges, "Charge")

# In[23]:


EventDisplayHist(digitimes, "Time")

# In[24]:


EventDisplayHist(digitimes, "Time", [940, 1100])


# In[25]:


# Add a 2d histogram of charge versus time


# In[26]:


def ChargeTimeHist(times, charges, title='Event Charge versus Time', cutrange=[[-1, -1], [-1, -1]]):
    """
    Makes a 2d histogram of charge versus time.
    inputs:
    times is an np.array of times of PMT hits
    charges is an np.array of charges of PMT hits
    title is the title of the histogram
    cutrange has two ranges, one in x and one in y [ [tmin, tmax], [qmin,qmax] ]
    """
    fig = plt.figure(figsize=[12, 12])
    tmin = times.min()
    tmax = times.max()
    qmin = charges.min()
    qmax = charges.max()

    if cutrange[0][0] != cutrange[0][1]:
        tmin = cutrange[0][0]
        tmax = cutrange[0][1]
    if cutrange[1][0] != cutrange[1][1]:
        qmin = cutrange[1][0]
        qmax = cutrange[1][1]

    plt.hist2d(times, charges, [100, 100], [[tmin, tmax], [qmin, qmax]])
    fig.suptitle(title, fontsize=20)
    plt.xlabel('Time (ns)', fontsize=18)
    plt.ylabel('Charge (pe)', fontsize=16)
    # plt.set_cmap('gist_heat_r')
    plt.set_cmap('cubehelix_r')
    plt.colorbar()


# In[27]:


ChargeTimeHist(digitimes, digicharges)

# In[28]:


ChargeTimeHist(digitimes, digicharges, 'QT', [[940, 1040], [0, 20]])


# # Pick a random event to display

# In[36]:


def GetRandomEvent():
    evno = random.randint(0, number_of_events - 1)
    dgitubes = datafile['digi_hit_pmt'][evno]
    dgicharges = datafile['digi_hit_charge'][evno]
    dgitimes = datafile['digi_hit_time'][evno]
    return (evno, dgitubes, dgicharges, dgitimes)


# def GetEvent( evnum ):
#    dgitubes = datafile[ 'digi_hit_pmt' ][ evnum ]
#    dgicharges = datafile[ 'digi_hit_charge' ][ evnum ]
#    dgitimes = datafile[ 'digi_hit_time' ][ evnum ]    
#    return (evnum, dgitubes, dgicharges, dgitimes)


# Actually, first let's draw event 124 which has a PMT in the event with PMT number 0

# In[39]:


evno = 124
digitubes, digicharges, digitimes = GetEvent(evno)
print("Displaying event number ", evno)
EventDisplay(digitubes, digicharges, "Charges for event " + str(evno))

# In[31]:


evno, digitubes, digicharges, digitimes = GetRandomEvent()
print("Displaying event number ", evno)
EventDisplay(digitubes, digicharges, "Charges for event " + str(evno))


# # Try to display the PMTs hit in a 3D view

# In[32]:


def TubesToXYZ(geotubes, geoxyz, intubes):
    """ 
        First two inputs are the geometry np.arrays of tube numbers and corresponding xyz locations.
        
        Take an np.array of PMT tube numbers (intubes), and corresponding np.arrays of the x,y,z of those tubes
        in outx, outy and outz
        
    """
    outx = np.array([], dtype=np.int32)
    outy = np.array([], dtype=np.int32)
    outz = np.array([], dtype=np.int32)
    for tube in intubes:
        idx = tube - 1
        outx = np.append(outx, [geoxyz[idx, 0]])
        outy = np.append(outy, [geoxyz[idx, 1]])
        outz = np.append(outz, [geoxyz[idx, 2]])

    return (outx, outy, outz)


# In[33]:


def GetParticleStartStops(datain, evno):
    """
    Return start and stop of each non-zero charged particle type
    
    ( x,  y , z, pid, Energy)
    
    where x -> [ [startx, stopx ] , [startx, stopx ], ... ]
    ...

    """
    xret = []
    yret = []
    zret = []
    pidret = []
    energies = []

    startpos = datain['track_start_position'][evno]
    stoppos = datain['track_stop_position'][evno]
    pids = datain['track_pid'][evno]
    enes = datain['track_energy'][evno]

    for idx, pid in enumerate(pids):
        # keep e, mu, gamma (11,13,22)
        if idx == 0:
            continue
        if abs(pid) != 11 and abs(pid) != 13 and pid != 22:
            continue
        xret.append([startpos[idx, 0], stoppos[idx, 0]])
        yret.append([startpos[idx, 1], stoppos[idx, 1]])
        zret.append([startpos[idx, 2], stoppos[idx, 2]])
        pidret.append(pid)
        energies.append(enes[idx])

    return (xret, yret, zret, pidret, energies)


# In[34]:


gstartpos = datafile['track_start_position'][evno]
gstoppos = datafile['track_stop_position'][evno]
gpids = datafile['track_pid'][evno]
genes = datafile['track_energy'][evno]
print(gstartpos)
print(gstoppos)
print(gpids)
print(genes)
partx, party, partz, partid, partene = GetParticleStartStops(datafile, evno)
print(partx)
print(party)
print(partz)
print(partid)

# In[139]:


trigger_time = datafile['trigger_time']
print("trigger_time=", trigger_time)
for ev in range(0, 100):
    first_trigger = np.argmin(trigger_time[ev])
    print("ev=", ev, " ntriggers=", len(trigger_time[ev]), " first_trigger=", first_trigger)

# In[141]:


digitubes = datafile['digi_hit_pmt']
decaytimes = np.array([])
for trigger_times in trigger_time:
    if len(trigger_times) == 2:
        decaytimes = np.append(decaytimes, trigger_times[1] - trigger_times[0])

decaytimes
EventDisplayHist(decaytimes, "Decay time (ns)")

# In[57]:


digitubes0 = datafile['digi_hit_pmt'][0]
digitubes0

# In[54]:


evno = 8
partx, party, partz, partid, partene = GetParticleStartStops(datafile, evno)
digitubes, digicharges, digitimes = GetEvent(evno)
evx, evy, evz = TubesToXYZ(tubes, tube_xyz, digitubes)
# evxyz = tube_xyz[ digitubes ]

digi_hit_trigger = datafile['digi_hit_trigger']
trigger_time = datafile['trigger_time']
print("trigger_time=", trigger_time)
first_trigger = np.argmin(trigger_time[evno])
print("first_trigger=", first_trigger)

good_hits1 = np.where(digi_hit_trigger[evno] == first_trigger)
hit_pmts1 = digitubes[good_hits1]
hit_charges1 = digicharges[good_hits1]
evxyz1 = tube_xyz[hit_pmts1]

second_trigger = np.argmax(trigger_time[evno])
good_hits2 = np.where(digi_hit_trigger[evno] == second_trigger)
hit_pmts2 = digitubes[good_hits2]
hit_charges2 = digicharges[good_hits2]
evxyz2 = tube_xyz[hit_pmts2]

print("trigger time 1 = ", trigger_time[evno][first_trigger], " trigger time 2 = ", trigger_time[evno][second_trigger])
print("allhits = ", len(digitubes), " first trigger hits=", len(hit_pmts1), " second trigger hits=", len(hit_pmts2))

# hit_pmts = datafile['true_hit_pmt'][evno]
# hit_positions = geofile['position'][hit_pmts]

fig = plt.figure(figsize=[10, 8])
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d(-R, R)
ax.set_ylim3d(-H / 2, H / 2)
ax.set_zlim3d(-R, R)
# ax.scatter( hit_positions[:,0], hit_positions[:,1], hit_positions[:,2], marker='.', c='b' )
# img = ax.scatter( evx, evy, evz, marker='o', c=digitimes-945, s=digicharges, vmin=940 , vmax=1040  )
# img = ax.scatter( evx, evy, evz, marker='o', c=digicharges, s=3 )#, vmin=940 , vmax=1040  )

img = ax.scatter(evxyz1[:, 0], evxyz1[:, 1], evxyz1[:, 2], marker='o', c=hit_charges1, s=3)

ax.scatter(evxyz2[:, 0], evxyz2[:, 1], evxyz2[:, 2], marker='.', c='r', s=3)

# plt.set_cmap('hot_r')
plt.set_cmap('cubehelix_r')
colors = {11: 'r', 13: 'g', 22: 'c'}
for i, pid in enumerate(partid):
    ax.plot(partx[i], party[i], partz[i], c=colors[abs(pid)])
    ax.text(partx[i][0], party[i][0], partz[i][0], "%.1f MeV" % partene[i], color=colors[abs(pid)])
ax.set_xlabel('X (cm)')
ax.set_ylabel('Y (cm)')
ax.set_zlabel('Z (cm)')
ax.set_title('Event %d' % evno)
# ax.view_init(elev=45.0, azim=45.0)

# plt.set_cmap('hot_r')
plt.colorbar(img)
plt.show()

# rotate the axes and update
# for angle in range(0, 360):
#    ax.view_init(30, angle)
#    plt.draw()
#    plt.pause(.001)


# In[ ]:


# In[ ]:


# In[84]:


import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure(figsize=[10, 8])
ax = fig.add_subplot(111, projection='3d')

f = np.load("event998.npz", allow_pickle=True)
g = np.load("mpmt_full_geo.npz")
truehit_pmts = f['true_hit_pmt'][1]
truehit_positions = g['position'][truehit_pmts]

truehit_parent = f['true_hit_parent'][1]
noisehit_pmts = truehit_pmts[truehit_parent == -1]
signalhit_pmts = truehit_pmts[truehit_parent != -1]
print(noisehit_pmts)
noisehit_positions = g['position'][noisehit_pmts]
signalhit_positions = g['position'][signalhit_pmts]
print(noisehit_positions)

hit_pmts = f['digi_hit_pmt'][1]
hit_positions = g['position'][hit_pmts]
hit_charges = f['digi_hit_charge'][1]

s_hit_pmts = np.array([], dtype=np.int32)
s_hit_charges = np.array([])
for idx, pmt in enumerate(hit_pmts):
    if np.isin(pmt, truehit_pmts):
        s_hit_pmts = np.append(s_hit_pmts, pmt)
        s_hit_charges = np.append(s_hit_charges, hit_charges[idx])

s_hit_positions = g['position'][s_hit_pmts]

ax = fig.add_subplot(111, projection='3d')
# ax.scatter(truehit_positions[:,0], truehit_positions[:,1], truehit_positions[:,2], marker='.')

# ax.scatter(hit_positions[:,0], hit_positions[:,1], hit_positions[:,2], marker=".")

# ax.scatter(noisehit_positions[:,0], noisehit_positions[:,1], noisehit_positions[:,2], marker="o", c='r')

ax.scatter(s_hit_positions[:, 0], s_hit_positions[:, 1], s_hit_positions[:, 2], marker="o", c=s_hit_charges, vmax=10)

plt.show()


# # Make new dataset that has random number of muons in each event
# 
# This dataset will be stored in an h5 file to match the format used by the convolutional neural network codes.  Existing single muon copy of npz file to h5 file have the format:
# 
# 
#     dtype_events = np.dtype(np.float32)
#     dtype_labels = np.dtype(np.int32)
#     dtype_energies = np.dtype(np.float32)
#     dtype_positions = np.dtype(np.float32)
#     dtype_IDX = np.dtype(np.int32)
#     dtype_PATHS = h5py.special_dtype(vlen=str)
#     dtype_angles = np.dtype(np.float32)
#     
#     h5_file = h5py.File(config.output_file[0], 'w')
#     
#     dset_event_data = h5_file.create_dataset("event_data",
#                                        shape=(num_nonzero_events,)+IMAGE_SHAPE,
#                                        dtype=dtype_events)
#                                        
#     dset_labels = h5_file.create_dataset("labels",
#                                    shape=(num_nonzero_events,),
#                                    dtype=dtype_labels)
#                                    
#     dset_energies = h5_file.create_dataset("energies",
#                                      shape=(num_nonzero_events, 1),
#                                      dtype=dtype_energies)
#     dset_positions = h5_file.create_dataset("positions",
#                                       shape=(num_nonzero_events, 1, 3),
#                                       dtype=dtype_positions)
#     dset_IDX = h5_file.create_dataset("event_ids",
#                                 shape=(num_nonzero_events,),
#                                 dtype=dtype_IDX)
#     dset_PATHS = h5_file.create_dataset("root_files",
#                                   shape=(num_nonzero_events,),
#                                   dtype=dtype_PATHS)
#     dset_angles = h5_file.create_dataset("angles",
#                                  shape=(num_nonzero_events, 2),
#                                  dtype=dtype_angles)
# 
# 
# * 'event_id' -> new event id, just increment a counter
# * 'root_file' -> use first muon's root file
# * 'pid'       -> ???
#  'position',
#  'direction',
#  'energy',
#  'digi_hit_pmt',
#  'digi_hit_charge',
#  'digi_hit_time',
#  'digi_hit_trigger',
#  'true_hit_pmt',
#  'true_hit_time',
#  'true_hit_pos',
#  'true_hit_start_time',
#  'true_hit_start_pos',
#  'true_hit_parent',
#  'track_id',
#  'track_pid',
#  'track_start_time',
#  'track_energy',
#  'track_start_position',
#  'track_stop_position',
#  'track_parent',
#  'track_flag',
#  'trigger_time']

# In[86]:


def GetNoise_and_Signal_PMTLists(evnum, data):
    """
    Inputs: 
        evnum == index into data for the event to get indices for
        data == npz datafile
        geo == npz geometry file
    
    Return two-tuple np.arrays:  one with the noise hits, and one with the non-noise hits
    """

    truehit_pmts = data['true_hit_pmt'][evnum]
    truehit_parent = data['true_hit_parent'][evnum]
    noisehit_pmts = truehit_pmts[truehit_parent == -1]
    signalhit_pmts = truehit_pmts[truehit_parent != -1]
    return (noisehit_pmts, signalhit_pmts)


def PMTChargeTime_in_list(evnum, pmtlist, data):
    """
    Return a tuple containing:
        np.array of pmt-numbers from evnum in data, that is in list of pmts
        np.array of digi-charges from evnum in data, that is in list of pmts
        np.array of digi-times from evnum in data, that is in list of pmts
    """
    allpmts = data['digi_hit_pmt'][evnum]
    allcharges = data['digi_hit_charge'][evnum]
    alltimes = data['digi_hit_time'][evnum]

    s_hit_pmts = np.array([], dtype=np.int32)
    s_hit_charges = np.array([])
    s_hit_times = np.array([])

    for idx, pmt in enumerate(allpmts):
        if np.isin(pmt, pmtlist):
            s_hit_pmts = np.append(s_hit_pmts, pmt)
            s_hit_charges = np.append(s_hit_charges, allcharges[idx])
            s_hit_times = np.append(s_hit_times, alltimes[idx])

    return (s_hit_pmts, s_hit_charges, s_hit_times)


# In[98]:


def SumEvents(event_numbers, time_offsets, datafile, only_noise=False):
    """
    This function sums the events in the list of event indices into the datafile.
    
    Inputs:
    event_numbers is list of indices into events in datfile to sum be summed
    time_offsets is list of time offests for each event (only used for >1 event)
    if only_noise true, then only return the noise hits from the first event in the list.
    
    Return tuple of ( tube, charge, time)
    tube is np.array of tubes (unique without duplicates)
    charge is np.array of charges (summed from all events hitting each corresponding tube)
    time is np.array of times (earliest from all events hitting each corresponding tube)
    
    Note: only take the noise hits from the first event in the list of event numbers.  
    """

    tube = np.array([], dtype=np.int32)
    charge = np.array([])
    time = np.array([])

    nev = len(event_numbers)
    if nev == 0:
        return (tube, charge, time)

        # Start with first event in list
    ievnum = event_numbers[0]
    if only_noise:
        noisepmts, signalpmts = GetNoise_and_Signal_PMTLists(ievnum, datafile)
        tube, charge, time = PMTChargeTime_in_list(ievnum, noisepmts, datafile)
    else:
        tube = datafile['digi_hit_pmt'][ievnum]
        charge = datafile['digi_hit_charge'][ievnum]
        time = datafile['digi_hit_time'][ievnum]

        # For remaining events only look at signal events
    for iev in range(1, len(event_numbers)):
        ievnum = event_numbers[iev]
        toffset = time_offsets[iev]
        noisepmts, signalpmts = GetNoise_and_Signal_PMTLists(ievnum, datafile)
        curtube, curcharge, curtime = PMTChargeTime_in_list(ievnum, signalpmts, datafile)
        for idx, pmt in enumerate(curtube):
            if np.isin(pmt, tube):
                for jdx, jpmt in enumerate(tube):
                    if jpmt == pmt:
                        charge[jdx] += curcharge[idx]
                        if curtime[idx] + toffset < time[jdx]:
                            time[jdx] = curtime[idx] + toffset
            else:
                # new pmt not in list yet... append
                tube = np.append(tube, pmt)
                charge = np.append(charge, curcharge[idx])
                time = np.append(time, curtime[idx] + toffset)

    return (tube, charge, time)


# In[109]:


# hit_charges1 = digicharges[good_hits]
# evxyz1 = tube_xyz[ hit_pmts1 ]


evlist = [50, 100, 150]

sumtubes, sumcharges, sumtimes = SumEvents(evlist, [0., 0., 0.], datafile)

print("total energy = ", np.sum(datafile['energy'][evlist]))
total_energy = 0
energies = datafile['energy']
for ev in evlist:
    print("event ", ev, " energy is ", energies[ev])
    total_energy += energies[ev]
print("energy=", total_energy)

evxyz = tube_xyz[sumtubes]

fig = plt.figure(figsize=[10, 8])
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d(-R, R)
ax.set_ylim3d(-H / 2, H / 2)
ax.set_zlim3d(-R, R)
# ax.scatter( hit_positions[:,0], hit_positions[:,1], hit_positions[:,2], marker='.', c='b' )
# img = ax.scatter( evx, evy, evz, marker='o', c=digitimes-945, s=digicharges, vmin=940 , vmax=1040  )
# img = ax.scatter( evx, evy, evz, marker='o', c=digicharges, s=3 )#, vmin=940 , vmax=1040  )

img = ax.scatter(evxyz[:, 0], evxyz[:, 1], evxyz[:, 2], marker='o', c=sumtimes, s=sumcharges, vmin=940, vmax=1040)

ax.set_xlabel('X (cm)')
ax.set_ylabel('Y (cm)')
ax.set_zlabel('Z (cm)')
ax.set_title('Event %d' % evno)
# ax.view_init(elev=45.0, azim=45.0)

# plt.set_cmap('hot_r')
plt.colorbar(img)
plt.show()


# In[129]:
def foo(mylist):
    item0 = mylist[0]


datafile['root_file']
datafile['position'][0]

# In[136]:


import os
import h5py
import numpy as np

IMAGE_SHAPE = (40, 40, 38)
PMT_LABELS = "PMTlabelSheet3.csv"


def count_events(files):
    # Because we want to remove events with 0 hits, 
    # we need to count the events beforehand (to create the h5 file).
    # This function counts and indexes the events with more than 0 hits.
    # Files need to be iterated in the same order to use the indexes.
    """ This is where we manually specify the file"""
    num_events = 0
    nonzero_file_events = []
    for file_index, f in enumerate(files):
        data = np.load(f, allow_pickle=True)
        nonzero_file_events.append([])
        hits = data['digi_hit_pmt']
        for i in range(len(hits)):
            if len(hits[i]) != 0:
                nonzero_file_events[file_index].append(i)
                num_events += 1
    return num_events, nonzero_file_events


def GenMapping(csv_file):
    mPMT_to_index = {}
    with open(csv_file) as f:
        rows = f.readline().split(",")[1:]
        rows = [int(r.strip()) for r in rows]

        for line in f:
            line_split = line.split(",")
            col = int(line_split[0].strip())
            for row, value in zip(rows, line_split[1:]):
                value = value.strip()
                if value:  # If the value is not empty
                    mPMT_to_index[int(value)] = [col, row]
    npmap = np.zeros((max(mPMT_to_index) + 1, 2), dtype=np.int)
    for k, v in mPMT_to_index.items():
        npmap[k] = v
    return npmap


def GenerateMultiMuonSample_h5(avg_mu_per_ev=2.5, sigma_time_offset=21.2):
    """
    Inputs: 
     avg_mu_per_ev == Poisson distribution mean for number of muons in each spill
     sigma_time_offset == Width of spill (Gaussian) in nanoseconds
    """
    files = ['event998.npz']

    # Remove whitespace 
    files = [x.strip() for x in files]

    # Check that files were provided
    if len(files) == 0:
        raise ValueError("No files provided!!")
    print("Merging " + str(len(files)) + " files")

    # Start merging
    num_nonzero_events, nonzero_event_indexes = count_events(files)
    print(num_nonzero_events)

    # np.random.poisson( avg_mu_per_ev, number_of_throws )
    num_muons = np.random.poisson(avg_mu_per_ev, num_nonzero_events)

    #

    dtype_events = np.dtype(np.float32)
    dtype_labels = np.dtype(np.int32)
    dtype_energies = np.dtype(np.float32)
    dtype_positions = np.dtype(np.float32)
    dtype_IDX = np.dtype(np.int32)
    dtype_PATHS = h5py.special_dtype(vlen=str)
    dtype_angles = np.dtype(np.float32)
    h5_file = h5py.File('multimuonfile.h5', 'w')
    dset_event_data = h5_file.create_dataset("event_data",
                                             shape=(num_nonzero_events,) + IMAGE_SHAPE,
                                             dtype=dtype_events)
    dset_labels = h5_file.create_dataset("labels",
                                         shape=(num_nonzero_events,),
                                         dtype=dtype_labels)
    dset_energies = h5_file.create_dataset("energies",
                                           shape=(num_nonzero_events, 1),
                                           dtype=dtype_energies)
    dset_positions = h5_file.create_dataset("positions",
                                            shape=(num_nonzero_events, 1, 3),
                                            dtype=dtype_positions)
    dset_IDX = h5_file.create_dataset("event_ids",
                                      shape=(num_nonzero_events,),
                                      dtype=dtype_IDX)
    dset_PATHS = h5_file.create_dataset("root_files",
                                        shape=(num_nonzero_events,),
                                        dtype=dtype_PATHS)
    dset_angles = h5_file.create_dataset("angles",
                                         shape=(num_nonzero_events, 2),
                                         dtype=dtype_angles)

    # 22 -> gamma, 11 -> electron, 13 -> muon
    # corresponds to labelling used in CNN with only barrel
    # IWCDmPMT_4pi_full_tank_gamma_E0to1000MeV_unif-pos-R371-y521cm_4pi-dir_3000evts_329.npz has an event
    # with pid 11 though....
    # pid_to_label = {22:0, 11:1, 13:2}

    offset = 0
    offset_next = 0
    mPMT_to_index = GenMapping(PMT_LABELS)
    # Loop over files
    for file_index, filename in enumerate(files):
        data = np.load(filename, allow_pickle=True)
        nonzero_events_in_file = len(nonzero_event_indexes[file_index])
        x_data = np.zeros((nonzero_events_in_file,) + IMAGE_SHAPE,
                          dtype=dtype_events)
        digi_hit_pmt = data['digi_hit_pmt']
        # digi_hit_charge = data['digi_hit_charge']
        # digi_hit_time = data['digi_hit_time']
        # digi_hit_trigger = data['digi_hit_trigger']
        # trigger_time = data['trigger_time']
        delay = 0
        # Loop over events in file
        # Loop over number of muons in each event
        event_id = np.array([], dtype=np.int32)
        root_file = np.array([], dtype=np.str)
        pid = np.array([])
        position = np.array([])
        direction = np.array([])
        energy = np.array([])
        labels = np.array([])

        for i, nmu in enumerate(num_muons):
            print("processing output entry ", i, " with ", nmu, " muons")
            indices = np.random.randint(0, len(digi_hit_pmt), max(1, nmu))
            time_offs = [0.]
            if nmu > 1:
                time_offs = np.append(time_offs, np.random.normal(0., sigma_time_offset, nmu - 1))
            hit_pmts, charge, time = SumEvents(indices, time_offs, data, nmu == 0)
            hit_mpmts = hit_pmts // 19
            pmt_channels = hit_pmts % 19
            rows = mPMT_to_index[hit_mpmts, 0]
            cols = mPMT_to_index[hit_mpmts, 1]
            x_data[i - delay, rows, cols, pmt_channels] = charge
            x_data[i - delay, rows, cols, pmt_channels + 19] = time

            # fix below!!!
            idx0 = indices[0]
            event_id = np.append(event_id, data['event_id'][idx0])
            root_file = np.append(root_file, data['root_file'][idx0])
            pid = np.append(pid, data['pid'][idx0])
            position = np.append(position, data['position'][idx0])
            direction = np.append(direction, data['direction'][idx0])
            energy = np.append(energy, np.sum(data['energy'][indices]))
            labels = np.append(labels, nmu)

        offset_next += nonzero_events_in_file

        file_indices = nonzero_event_indexes[file_index]

        dset_IDX[offset:offset_next] = event_id[file_indices]
        dset_PATHS[offset:offset_next] = root_file[file_indices]
        dset_energies[offset:offset_next, :] = energy[file_indices].reshape(-1, 1)
        dset_positions[offset:offset_next, :, :] = position[file_indices].reshape(-1, 1, 3)
        dset_labels[offset:offset_next] = labels[file_indices]

        direction = direction[file_indices]
        polar = np.arccos(direction[:, 1])
        azimuth = np.arctan2(direction[:, 2], direction[:, 0])
        dset_angles[offset:offset_next, :] = np.hstack((polar.reshape(-1, 1), azimuth.reshape(-1, 1)))
        dset_event_data[offset:offset_next, :] = x_data

        offset = offset_next
        print("Finished file: {}".format(filename))

    print("Saving")
    h5_file.close()
    print("Finished")


# In[ ]:


GenerateMultiMuonSample_h5(avg_mu_per_ev=2.5, sigma_time_offset=21.2)
