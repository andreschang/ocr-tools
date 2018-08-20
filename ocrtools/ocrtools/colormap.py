# -*- coding: utf-8 -*-
## Andres Chang 
## 2017-12-28
## Code to make color maps from
## a list of hex colors
## 12-28 updated with cmap_sequence
## code to make a dynamic series of cmaps
##
## also split make_cmap into
## dedicated mapping function (make_cmap)
## that calls from dedicated make_grad function
## in 04.01, edited make_grad function to
## include not even spacing
##
## 12-30 added transition function that is
## automatically called if the step between
## two mains in adjacent cmaps switch directions
##
## 12-31 cleaned up
## updated to lower saturation during transition
##
## 01-01
## updated with ngl (non-object) cmap output
##
## 08-11
## updated to specify data-types and print() for python 3

import numpy as np
import matplotlib.pyplot as plt
from colour import *
from matplotlib.colors import LinearSegmentedColormap

class colormap(object):

  def __init__(self):
    self.name = 'colormap'
    

  def make_grad(colors, nsteps = 400, spacing = 'even', output = 'hex', rv_fix = 'on',\
   name = 'custom_grad', sequence = 'off', ngl = 'off'):
    ## colors is a list of hex values gradient is filled in between
    ## default is 400 steps with even spacing
    ## but spacing can be set to a list of main color indices
    ## for uneven spacing of main colors

    nmains = int(len(colors))
    nspans = nmains-1 # number of grad spans
    n = int(nsteps/nspans)
    nsteps = int(nsteps)

    mains = []
    hsl_mains = np.zeros((nmains, 3))
    hsl_array = np.zeros((nsteps, 3))

    for i in range(nmains):
      mains.append(Color(colors[i]))
      hsl_mains[i] = [mains[i].hue, mains[i].saturation, mains[i].luminance]

    hsl_step = np.diff(hsl_mains, axis = 0)
    hsl_step = hsl_step.tolist()

    if spacing == 'even':
      for i in range(nspans):
        # if i == 1:
          # print name
          # print hsl_mains[i][0], hsl_mains[i+1][0]
          # print hsl_step[i][0]
        for j in range(3):
           ## for hue transitions that cross the between red and violet
          if (rv_fix == 'on' and j == 0 and np.abs(hsl_step[i][j]) > 0.5):
            if hsl_step[i][j] > 0.5:
              new_step = -(1.-hsl_step[i][j])
              hsl_array[i*n:(i+1)*n, j] = np.concatenate((np.arange(hsl_mains[i][j],0.0, new_step/n), \
                np.arange(1.0, hsl_mains[i+1][j], new_step/n)))[:n]
            elif hsl_step[i][j] < -0.5:
              new_step = -(-1.-hsl_step[i][j])
              hsl_array[i*n:(i+1)*n, j] = np.concatenate((np.arange(hsl_mains[i][j],1.0, new_step/n), \
                np.arange(0.0, hsl_mains[i+1][j], new_step/n)))[:n]
            ## -- for hue transitions that cross the between red and violet
          else:
            if hsl_step[i][j] == 0:
              hsl_array[i*n:(i+1)*n, j] = hsl_mains[i][j]*np.ones(n)
            else:
              hsl_array[i*n:(i+1)*n, j] = np.arange(hsl_mains[i][j], hsl_mains[i+1][j], hsl_step[i][j]/n)[:n]

    else:
      n = np.diff(spacing)
      ni = 0
      for i in range(nspans):
        for j in range(3):
          ## for hue transitions that cross the between red and violet
          if (rv_fix == 'on' and j == 0 and np.abs(hsl_step[i][j]) > 0.5):
            if hsl_step[i][j] > 0.5:
              new_step = -(1.-hsl_step[i][j])
              hsl_array[ni:ni+n[i], j] = np.concatenate((np.arange(hsl_mains[i][j],0.0, new_step/n[i]), \
                np.arange(1.0, hsl_mains[i+1][j], new_step/n[i])))[:n[i]]
            elif hsl_step[i][j] < -0.5:
              new_step = -(-1.-hsl_step[i][j])
              hsl_array[ni:ni+n[i], j] = np.concatenate((np.arange(hsl_mains[i][j],1.0, new_step/n[i]), \
                np.arange(0.0, hsl_mains[i+1][j], new_step/n[i])))[:n[i]]
            ## -- for hue transitions that cross the between red and violet
          else:
            hsl_array[ni:ni+n[i], j] = np.arange(hsl_mains[i][j], hsl_mains[i+1][j], hsl_step[i][j]/n[i])[:n[i]]

        ni += n[i]

    if output == 'hex':
      hsl_array = np.array([Color( hsl = (hsl_array[nn, 0], hsl_array[nn, 1], hsl_array[nn, 2])).hex for \
      nn in range(nsteps)])
    
    elif output == 'Color':
      hsl_array = np.array([Color( hsl = (hsl_array[nn, 0], hsl_array[nn, 1], hsl_array[nn, 2])) for \
      nn in range(nsteps)])

    ## if output == 'hsl':
    ## do nothing to the array

    if sequence == 'off':
      return hsl_array
    elif sequence == 'on':
      return hsl_array, np.array(hsl_step)[:,0], hsl_mains


  def make_cmap(colors, nsteps = 400, name = 'custom_map', rv_fix = 'on', edit = '', sequence = 'off', ngl = 'off'):

    if ngl == "on":

      hex_array = colormap.make_grad(colors, nsteps = nsteps, output = 'hex', rv_fix = rv_fix, name = name, sequence = 'off')

      if sequence == 'on':
        hsl_array, step1, hsls = colormap.make_grad(colors, nsteps = nsteps, output = 'hsl', rv_fix = rv_fix, name = name, sequence = 'on')

      newmap = np.zeros((nsteps, 3))
      for c in range(nsteps):
        rgb = list(hex2rgb(hex_array[c]))
        newmap[c] = rgb

    else:

      if sequence == "off":
        hsl_array = colormap.make_grad(colors, nsteps = nsteps, output = 'hsl', rv_fix = rv_fix, name = name, sequence = 'off')
      else:
        hsl_array, step1, hsls = colormap.make_grad(colors, nsteps = nsteps, output = 'hsl', rv_fix = rv_fix, name = name, sequence = 'on')

      cdict = {'red':[], 'green':[], 'blue':[]}

      for i in range(nsteps):
        cpick = Color(hsl = (hsl_array[i, 0], hsl_array[i, 1], hsl_array[i, 2]))
        cdict['red'].append((float(i)/nsteps, cpick.red, cpick.red))
        cdict['green'].append((float(i)/nsteps, cpick.green, cpick.green))
        cdict['blue'].append((float(i)/nsteps, cpick.blue, cpick.blue))

      cdict['red'].append((1., cdict['red'][-1][1], cdict['red'][-1][2]))
      cdict['green'].append((1., cdict['green'][-1][1], cdict['green'][-1][2]))
      cdict['blue'].append((1., cdict['blue'][-1][1], cdict['blue'][-1][2]))

      cdict['red'], cdict['green'], cdict['blue'] = tuple(cdict['red']), tuple(cdict['green']), tuple(cdict['blue'])
      newmap = LinearSegmentedColormap(name, cdict)

    if sequence == "off":
      return newmap
    else:
      return newmap, step1, hsls

  ## color_lists is a list of each combination of colors (high, mid, low)
  ## ex. [[winter colors], [spring colors], [summer colors], [spring colors]]

  def cmap_sequence(all_mains, spacing = 'even', nmaps = 365, cycle = 'yes', rv_fix1 = "on", rv_fix2 = "off", ngl = 'off'):

    ## all_mains is a list of each combination of colors (high, mid, low)
    ## ex. [[winter colors], [spring colors], [summer colors], [spring colors]]

    ## add first color to end of list
    ## so that cmap is cyclical
    if cycle == 'yes':
      all_mains.append(all_mains[0])

    ngroups = len(all_mains)
    nmaps = int(nmaps)

    ## check to see if each main list has same num of mains
    ## if not, calculate gradient for those with fewer and
    ## all colors in between first and last mains are redefined
    ## to match number and spacing of num_mains

    num_mains = np.max([len(sublist) for sublist in all_mains])

    for sublist in range(ngroups):
      if len(all_mains[sublist]) < num_mains:
        ninsert = num_mains - 2
        fullgrad = colormap.make_grad(all_mains[sublist], output = 'hsl', rv_fix = 'off')
        midpoints = [(n+1)*int(fullgrad.shape[0]/(ninsert+1)) for n in range(ninsert)]
        print('ADDING MAINS at ... '+str(midpoints))
        mcolors = [Color(hsl = (fullgrad[midpoints[nn], 0], fullgrad[midpoints[nn], 1],\
         fullgrad[midpoints[nn], 2])) for nn in range(ninsert)]
        all_mains[sublist] = [all_mains[sublist][0], all_mains[sublist][-1]]
        for n in range(ninsert):
          all_mains[sublist].insert((n+1), mcolors[n].hex)

    all_mains = np.array(all_mains)

    ## calculate horizontal gradients for main points
    ## (i.e. across time)
    ## then vertical gradients
    ## (i.e. each day)

    horizontals =  []
    for row in range(num_mains):
      row_mains = [cpick for cpick in all_mains[:, row]]
      horizontals.append(colormap.make_grad(row_mains, nsteps = nmaps, spacing = spacing, output = 'hex', rv_fix = rv_fix1))

    horizontals = np.swapaxes(np.array(horizontals), 0, 1)

    allmaps = []
    allhsls = []
    apply_when = -1

    for day in range(nmaps):
      if rv_fix2 == 'on':
        cmap, step1, hsls = colormap.make_cmap(horizontals[day], name = str(day)+'_map', rv_fix = 'on', sequence = 'on', ngl = ngl)
      else:
        cmap, step1, hsls = colormap.make_cmap(horizontals[day], name = str(day)+'_map', rv_fix = 'off', sequence = 'on', ngl = ngl)

      allhsls.append(hsls)

      if day == 0:
        step0 = step1

      for i in range(len(step1)):
        # print step0[i], step1[i]
        if step0[i]*step1[i] < -0.05:
          print("APPLY TRANSITION across break at "+str(day-1)+" - "+str(day))

          ## record which part of the gradient has reversal
          transitionl = 16
          rspan = i
          apply_when = day+transitionl/2

      if day == apply_when:
        sub_sequence = colormap.cmap_transition(allhsls[day-transitionl:day], rspan, transitionl, ngl = ngl)
        allmaps = allmaps[:day-transitionl]
        allmaps = allmaps + sub_sequence

      allmaps.append(cmap)

      step0 = step1

    return allmaps

  def cmap_transition(hsls, rspan, transitionl, ngl = 'off'):
    keyframe0 = int(transitionl/2-1)
    keyframe1 = int(transitionl/2)

    hsls[keyframe0][rspan] = hsls[keyframe0][rspan+1]
    hsls[keyframe0][rspan][1] = .7*hsls[keyframe0][rspan][1]
    hsls[keyframe0][rspan+1][1] = .7*hsls[keyframe0][rspan+1][1]
    hsls[keyframe1][rspan] = hsls[keyframe1][rspan+1]
    hsls[keyframe1][rspan][1] = .7*hsls[keyframe1][rspan][1]
    hsls[keyframe1][rspan+1][1] = .7*hsls[keyframe1][rspan+1][1]

    mainsb = [Color(hsl = (hsls[0][m][0], hsls[0][m][1], hsls[0][m][2])).hex for m in range(len(hsls[0]))]
    mains0 = [Color(hsl = (hsls[keyframe0][m][0], hsls[keyframe0][m][1], \
      hsls[keyframe0][m][2])).hex for m in range(len(hsls[0]))]
    mains1 = [Color(hsl = (hsls[keyframe1][m][0], hsls[keyframe1][m][1], \
      hsls[keyframe1][m][2])).hex for m in range(len(hsls[0]))]
    mainsf = [Color(hsl = (hsls[-1][m][0], hsls[-1][m][1], hsls[-1][m][2])).hex for m in range(len(hsls[0]))]
    all_mains0 = np.concatenate(([mainsb], [mains0]))
    all_mains1 = np.concatenate(([mains1], [mainsf]))

    trans0 = colormap.cmap_sequence(all_mains0, nmaps = transitionl/2, cycle = 'no', rv_fix2 = 'off', ngl = ngl)
    trans1 = colormap.cmap_sequence(all_mains1, nmaps = transitionl/2, cycle = 'no', ngl = ngl)
    trans_maps = trans0+trans1


    return trans_maps

  def test_cmap(cmap, fname = 'test_cmap'):
    x = np.arange(0,50,0.01)
    y = np.arange(0,50,0.01)
    xx, yy = np.meshgrid(x,y)
    Z = (xx**2+yy**2)
    fig, ax = plt.subplots(figsize = (3,3))
    ax.imshow(Z, cmap = cmap)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(fname+'.png', dpi = 160, bbox_inches='tight')
    plt.close()

