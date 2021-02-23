import matplotlib
import matplotlib.pyplot as plt
import numpy as np


labels = ['BC', 'POM', 'SO₄²⁻','SOA']
ENE = [-0.002857782716353798, -0.002038637294025899, -0.35557909605172633, -0.02135823121782214]
IND = [-0.15924662053629474, -0.2719267815157253, -0.8729732581946081, -0.08659916752816331]
RCO = [-0.12047899285243371, -0.5145862673284857, -0.0686171459935433, -0.021762166776834464]
SHP = [-0.011777480316432831, -0.037962471346245294, -0.0282736317512186, -0.023300925919002447]
SLV = [0.009041348385361383, 0.03065128070241038, 0.06256756541812161, 0.032692956784264475]
TRA = [-0.0008856104636625362, -0.020351704966094775, -0.015080816859415185, -0.02054822208624702]  
WST = [0.016403098092658847, 0.08156546728813334, 0.03389368558507823, 0.01604341018819297]

ENE_r = ENE/np.sum(ENE)
IND_r = IND/np.sum(IND)
RCO_r = RCO/np.sum(RCO)
SHP_r = SHP/np.sum(SHP)
SLV_r = SLV/np.sum(SLV)
TRA_r = TRA/np.sum(TRA)
WST_r = WST/np.sum(WST)

x = np.arange(len(labels))  # the label locations
width = 0.12  # the width of the bars
fig = plt.figure(figsize = (15,6))
ax = fig.add_subplot(1, 2, 1)
#fig, ax = plt.subplots()
rects1 = ax.bar(x-3*width, ENE_r, width, label='ENE')
rects2 = ax.bar(x-2*width, IND_r, width, label='IND')
rects3 = ax.bar(x-1*width, RCO_r, width, label='RCO')
rects4 = ax.bar(x, SHP_r, width, label='SHP')
rects5 = ax.bar(x+width, SLV_r, width, label='SLV')
rects6 = ax.bar(x+2*width, TRA_r, width, label='TRA')
rects7 = ax.bar(x+3*width, WST_r, width, label='WST')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('μg year⁻¹')
ax.set_title('Aerosol concentration variation in CN classified by sectors')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
# =============================================================================
#         ax.annotate('{}'.format(height),
#                     xy=(rect.get_x(), height),
#                     #xytext=(0, 3),  # 3 points vertical offset
#                     textcoords="offset points",
#                     fontsize = 8)
# =============================================================================
        if height>0:
            ax.text(rect.get_x(), height+0.2, '{}'.format(height))
        else:
            ax.text(rect.get_x(), height-1.2, '{}'.format(height))


autolabel(rects1)
autolabel(rects2)
autolabel(rects3)
autolabel(rects4)
autolabel(rects5)
autolabel(rects6)
autolabel(rects7)

plt.savefig('E:/test.jpeg',bbox_inches = 'tight')
plt.show()

#############################################################################
#
# ------------
#
# References
# """"""""""
#
# The use of the following functions, methods and classes is shown
# in this example:

matplotlib.axes.Axes.bar
matplotlib.pyplot.bar
matplotlib.axes.Axes.annotate
matplotlib.pyplot.annotate
