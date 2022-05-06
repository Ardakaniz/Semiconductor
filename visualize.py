import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter

def fn(arr):
    l = len(arr)
    new = np.array([arr for i in range(l)])
    return new

field = np.loadtxt(r"VS\semiconductor\electric_field.txt", dtype = float)
X = np.linspace(0, 1, 500);
print(field.shape)

data = np.loadtxt(r"VS\semiconductor\data.dat", dtype=float)
data = np.array(list(map(fn, data)))
#data = data[::10]

plt.rcParams["figure.autolayout"] = True

fig = plt.figure(dpi=120)
ax = fig.add_subplot(211)
field_ax = fig.add_subplot(212)
# I like to position my colorbars this way, but you don't have to
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')

lines,  = field_ax.plot([], [], c='r')
field_ax.set_xlim(0, 1)
field_ax.set_ylim(0, np.max(field) * 1.1)
field_ax.set_title("Electric field")
field_ax.set_xlabel("x")
field_ax.set_ylabel("E")

im = ax.imshow(data[0], cmap='seismic',interpolation="nearest",aspect="auto", extent=(0, 1, 0, 1))
ax.set_xlabel("x")
ax.set_title("Charge density")

colorbar_format = lambda x, pos: '+' if x > 0 else ('-' if x < 0 else '0')
cb = fig.colorbar(im,cax=cax,ax=ax,label='Charge density', format=FuncFormatter(colorbar_format))

def init():
    lines.set_data(X, field[0])
    return im, lines,

def animate(i):
    #cax.cla()
    im.set_data(data[i])
    field_ax.set_ylim(np.min(field[i]), np.max(field[i])*1.1)
    lines.set_ydata(field[i])
    #print(X[i], field[i,0])
    im.set_clim(np.min(data[i]), np.max(data[i]))
    #plt.colorbar(im, cax=cax)
    return im, lines,

plt.pause(5)
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(data), interval=1.0/30.0, blit=True)
plt.show()