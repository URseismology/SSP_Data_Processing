import matplotlib.pyplot as plt
from numpy import sin, cos, exp, pi, arange, mean, array
from matplotlib.transforms import Affine2D

# load the data
t = arange(0.0, 2.0, 0.01)
s1 = sin(2*pi*t)
s2 = exp(-t)
s3 = sin(2*pi*t)*exp(-t)
s4 = sin(2*pi*t)*cos(4*pi*t)
s5 = s1*s2
s6 = s1-s4
s7 = s3*s4-s1

signals = (s1, s2, s3, s4, s5, s6, s7)
for sig in signals:
    sig -= mean(sig)

lineprops = dict(linewidth=1, color='black', linestyle='-')
fig = plt.figure()
ax = fig.gca()

# scale -10 to 10 on the y scale.  -10 to 10 means that a signal with a min/max
# amplitude of 10 will span the entire vertical extent of the axes
scale = 10

ax.set_ylim(0, scale)

# now add the signals, set the transform, and set the offset of each line
ticklocs = []
for i, s in enumerate(signals):
    offset = (i+1.)/(len(signals)+1.)

    # The default transform is ax.transScale + (ax.transLimits + ax.transAxes),
    # but we are going to insert a slight translation in between the transLimits
    # step and the transAxes step.
    offsetTrans = ax.transScale + (ax.transLimits + Affine2D().translate(0, offset)
                                    + ax.transAxes)
    lines = ax.plot(t, s, transform=offsetTrans, **lineprops)
    ticklocs.append(offset)


ax.set_yticks(ticklocs)
ax.set_yticklabels(['S%d'%(i+1) for i in range(len(signals))])

# place all the y tick attributes in axes coords  
allticks = []
labels = []
ax.set_yticks(ticklocs)
for tick in ax.yaxis.get_major_ticks():
    allticks.extend(( tick.label1, tick.label2, tick.tick1line,
                      tick.tick2line, tick.gridline))
    labels.append(tick.label1)

plt.setp(allticks, transform=ax.transAxes)
plt.setp(labels, x=-0.01)

ax.set_xlabel('time (s)')


# Because we have hacked the transforms, you need a special method to
# set the voltage gain; this is a naive implementation of how you
# might want to do this in real life (eg make the scale changes
# exponential rather than linear) but it gives you the idea
def set_ygain(direction):
    set_ygain.scale += direction
    if set_ygain.scale <=0:
        set_ygain.scale -= direction
        return

    ax.set_ylim(0, set_ygain.scale)
    plt.draw()
set_ygain.scale = scale

def keypress(event):
    if event.key in ('+', '='): set_ygain(-1)
    elif event.key in ('-', '_'): set_ygain(1)

plt.connect('key_press_event', keypress)
ax.set_title('Use + / - to change y gain')
plt.show()

