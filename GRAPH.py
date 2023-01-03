import matplotlib.pyplot as plt
import matplotlib as mpl

# DOWNWELLING IRRADIANCE
# label axis
plt.title("Downwelling Irradiance", fontsize = 20, fontweight = "bold")
plt.xlabel("Wavelength (nm)", fontsize = 15)
# figure out how to program sub and superscripts
plt.ylabel("E", fontsize = 15)
mpl.rcParams["lines.linewidth"] = 3

# set range for axis
plt.xlim(400, 700)
plt.ylim(0.00, 0.14)

# create data arrays
x = [400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700]

y1 = [0.18, 0.24, 0.29, 0.32, 0.36, 0.39, 0.41, 0.43, 0.44, 0.44, 0.45, 0.41, 0.49]
y2 = [0.14, 0.18, 0.23, 0.26, 0.31, 0.35, 0.37, 0.39, 0.40, 0.39, 0.40, 0.35, 0.42]
y3 = [0.09, 0.13, 0.18, 0.21, 0.25, 0.30, 0.32, 0.33, 0.33, 0.33, 0.33, 0.29, 0.35]
y4 = [0.07, 0.10, 0.14, 0.17, 0.21, 0.25, 0.28, 0.30, 0.30, 0.29, 0.28, 0.24, 0.31]
y5 = [0.06, 0.08, 0.11, 0.14, 0.18, 0.22, 0.26, 0.27, 0.26, 0.26, 0.25, 0.20, 0.26]
y6 = [0.05, 0.07, 0.09, 0.12, 0.15, 0.20, 0.23, 0.25, 0.24, 0.24, 0.23, 0.17, 0.23]
y7 = [0.04, 0.05, 0.08, 0.10, 0.13, 0.17, 0.21, 0.23, 0.22, 0.22, 0.21, 0.15, 0.21]


# plot data on graph
plt.plot(x, y1, label = "0.00", color = "firebrick")
plt.plot(x, y2, label = "0.05", color = "darkorange")
plt.plot(x, y3, label = "0.10", color = "gold")
plt.plot(x, y4, label = "0.15", color = "mediumseagreen")
plt.plot(x, y5, label = "0.20", color = "mediumturquoise")
plt.plot(x, y6, label = "0.25", color = "cornflowerblue")
plt.plot(x, y7, label = "0.29", color = "mediumpurple")
plt.legend()
plt.show()


# UPWELLING IRRADIANCE
# label axis
plt.title("Upwelling Irradiance", fontsize = 20, fontweight = "bold")
plt.xlabel("Wavelength (nm)", fontsize = 15)
# figure out how to program sub and superscripts
plt.ylabel("E", fontsize = 15)
mpl.rcParams["lines.linewidth"] = 3

# set range for axis
plt.xlim(400, 700)
plt.ylim(0.00, 0.14)

# create data arrays
x = [400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700]

y1 = [0.18, 0.24, 0.29, 0.32, 0.36, 0.39, 0.41, 0.43, 0.44, 0.44, 0.45, 0.41, 0.49]
y2 = [0.14, 0.18, 0.23, 0.26, 0.31, 0.35, 0.37, 0.39, 0.40, 0.39, 0.40, 0.35, 0.42]
y3 = [0.09, 0.13, 0.18, 0.21, 0.25, 0.30, 0.32, 0.33, 0.33, 0.33, 0.33, 0.29, 0.35]
y4 = [0.07, 0.10, 0.14, 0.17, 0.21, 0.25, 0.28, 0.30, 0.30, 0.29, 0.28, 0.24, 0.31]
y5 = [0.06, 0.08, 0.11, 0.14, 0.18, 0.22, 0.26, 0.27, 0.26, 0.26, 0.25, 0.20, 0.26]
y6 = [0.05, 0.07, 0.09, 0.12, 0.15, 0.20, 0.23, 0.25, 0.24, 0.24, 0.23, 0.17, 0.23]
y7 = [0.04, 0.05, 0.08, 0.10, 0.13, 0.17, 0.21, 0.23, 0.22, 0.22, 0.21, 0.15, 0.21]


# plot data on graph
plt.plot(x, y1, label = "0.00", color = "firebrick")
plt.plot(x, y2, label = "0.05", color = "darkorange")
plt.plot(x, y3, label = "0.10", color = "gold")
plt.plot(x, y4, label = "0.15", color = "mediumseagreen")
plt.plot(x, y5, label = "0.20", color = "mediumturquoise")
plt.plot(x, y6, label = "0.25", color = "cornflowerblue")
plt.plot(x, y7, label = "0.29", color = "mediumpurple")
plt.legend()
plt.show()


# CANOPY REFLECTANCE
# label axis
plt.title("Canopy Reflectance", fontsize = 20, fontweight = "bold")
plt.xlabel("Wavelength (nm)", fontsize = 15)
plt.ylabel("Reflectance (rel)", fontsize = 15)
mpl.rcParams["lines.linewidth"] = 3

# set range for axis
plt.xlim(400, 700)
plt.ylim(0.00, 0.60)

# create data arrays
x = [400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700]

y1 = [0.18, 0.24, 0.29, 0.32, 0.36, 0.39, 0.41, 0.43, 0.44, 0.44, 0.45, 0.41, 0.49]
y2 = [0.14, 0.18, 0.23, 0.26, 0.31, 0.35, 0.37, 0.39, 0.40, 0.39, 0.40, 0.35, 0.42]
y3 = [0.09, 0.13, 0.18, 0.21, 0.25, 0.30, 0.32, 0.33, 0.33, 0.33, 0.33, 0.29, 0.35]
y4 = [0.07, 0.10, 0.14, 0.17, 0.21, 0.25, 0.28, 0.30, 0.30, 0.29, 0.28, 0.24, 0.31]
y5 = [0.06, 0.08, 0.11, 0.14, 0.18, 0.22, 0.26, 0.27, 0.26, 0.26, 0.25, 0.20, 0.26]
y6 = [0.05, 0.07, 0.09, 0.12, 0.15, 0.20, 0.23, 0.25, 0.24, 0.24, 0.23, 0.17, 0.23]
y7 = [0.04, 0.05, 0.08, 0.10, 0.13, 0.17, 0.21, 0.23, 0.22, 0.22, 0.21, 0.15, 0.21]


# plot data on graph
plt.plot(x, y1, label = "0.00", color = "firebrick")
plt.plot(x, y2, label = "0.05", color = "darkorange")
plt.plot(x, y3, label = "0.10", color = "gold")
plt.plot(x, y4, label = "0.15", color = "mediumseagreen")
plt.plot(x, y5, label = "0.20", color = "mediumturquoise")
plt.plot(x, y6, label = "0.25", color = "cornflowerblue")
plt.plot(x, y7, label = "0.30", color = "mediumpurple")
plt.legend()
plt.show()